// tAOFlaggerStep.cc: Test program for class AOFlaggerStep
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Ger van Diepen

#include "tStepCommon.h"
#include "mock/ThrowStep.h"

#include <boost/test/unit_test.hpp>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include "../../AOFlaggerStep.h"

#include <dp3/base/DP3.h>
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>

#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::steps::AOFlaggerStep;
using dp3::steps::Step;

using dp3::common::ParameterSet;

using casacore::Cube;
using casacore::MPosition;
using casacore::Quantum;

BOOST_AUTO_TEST_SUITE(aoflaggerstep)

// Simple class to generate input arrays.
// It can only set all flags to true or all to false.
// Weights are always 1.
// It can be used with different nr of times, channels, etc.
class TestInput : public dp3::steps::MockInput {
 public:
  TestInput(int ntime, int nant, int nchan, int ncorr, bool flag)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr),
        itsFlag(flag) {
    info() = DPInfo(itsNCorr, itsNChan);
    info().setTimes(100.0, 100.0 + (itsNTime - 1) * 5.0, 5.0);

    // Fill the baseline stations; use 4 stations.
    // So they are called 00 01 02 03 10 11 12 13 20, etc.
    std::vector<int> ant1(itsNBl);
    std::vector<int> ant2(itsNBl);
    int st1 = 0;
    int st2 = 0;
    for (int i = 0; i < itsNBl; ++i) {
      ant1[i] = st1;
      ant2[i] = st2;
      if (++st2 == 4) {
        st2 = 0;
        if (++st1 == 4) {
          st1 = 0;
        }
      }
    }
    std::vector<std::string> antNames(4);
    antNames[0] = "rs01.s01";
    antNames[1] = "rs02.s01";
    antNames[2] = "cs01.s01";
    antNames[3] = "cs01.s02";
    // Define their positions (more or less WSRT RT0-3).
    std::vector<casacore::MPosition> antPos(4);
    casacore::Vector<double> vals(3);
    vals[0] = 3828763;
    vals[1] = 442449;
    vals[2] = 5064923;
    antPos[0] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    vals[0] = 3828746;
    vals[1] = 442592;
    vals[2] = 5064924;
    antPos[1] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    vals[0] = 3828729;
    vals[1] = 442735;
    vals[2] = 5064925;
    antPos[2] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    vals[0] = 3828713;
    vals[1] = 442878;
    vals[2] = 5064926;
    antPos[3] = MPosition(Quantum<casacore::Vector<double>>(vals, "m"),
                          MPosition::ITRF);
    std::vector<double> antDiam(4, 70.);
    info().set(antNames, antDiam, antPos, ant1, ant2);

    // Define the frequencies.
    std::vector<double> chanFreqs(nchan);
    std::vector<double> chanWidth(nchan, 100000.);
    std::iota(chanFreqs.begin(), chanFreqs.end(), 1050000);
    info().setChannels(std::move(chanFreqs), std::move(chanWidth));
  }

 private:
  virtual bool process(const DPBuffer&) {
    // Stop when all times are done.
    if (itsCount == itsNTime) {
      return false;
    }
    Cube<std::complex<float>> data(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(data.size()); ++i) {
      data.data()[i] = std::complex<float>(1.6, 0.9);
    }
    if (itsCount == 5) {
      data += std::complex<float>(10., 10.);
    }
    DPBuffer buf;
    buf.setTime(itsCount * 5 + 2);  // same interval as in updateAveragInfo
    buf.setData(data);
    Cube<float> weights(data.shape());
    weights = 1.;
    buf.setWeights(weights);
    Cube<bool> flags(data.shape());
    flags = itsFlag;
    buf.setFlags(flags);
    // The fullRes flags are a copy of the XX flags, but differently shaped.
    // They are not averaged, thus only 1 time per row.
    Cube<bool> fullResFlags(itsNChan, 1, itsNBl);
    fullResFlags = itsFlag;
    buf.setFullResFlags(fullResFlags);
    getNextStep()->process(buf);
    ++itsCount;
    return true;
  }

  virtual void finish() { getNextStep()->finish(); }
  void updateInfo(const DPInfo&) override {
    // Do nothing / keep the info set in the constructor.
  }

  int itsCount, itsNTime, itsNBl, itsNChan, itsNCorr;
  bool itsFlag;
};

// Class to check result.
class TestOutput : public dp3::steps::test::ThrowStep {
 public:
  TestOutput(int ntime, int nant, int nchan, int ncorr)
      : itsCount(0),
        itsNTime(ntime),
        itsNBl(nant * (nant + 1) / 2),
        itsNChan(nchan),
        itsNCorr(ncorr) {}

 private:
  virtual bool process(const DPBuffer& buf) {
    // Fill expected result in similar way as TestInput.
    Cube<std::complex<float>> result(itsNCorr, itsNChan, itsNBl);
    for (int i = 0; i < int(result.size()); ++i) {
      result.data()[i] = std::complex<float>(1.6, 0.9);
    }
    if (itsCount == 5) {
      result += std::complex<float>(10., 10.);
    }
    // Check the result.
    /// cout << buf.getData()<< result;
    BOOST_CHECK(allNear(real(buf.getData()), real(result), 1e-10));
    BOOST_CHECK(allNear(imag(buf.getData()), imag(result), 1e-10));
    BOOST_CHECK_CLOSE(buf.getTime(), 2 + 5.0 * itsCount, 1e-8);
    ++itsCount;
    return true;
  }

  virtual void finish() {}
  virtual void updateInfo(const DPInfo& info) {
    BOOST_CHECK(int(info.origNChan()) == itsNChan);
    BOOST_CHECK(int(info.nchan()) == itsNChan);
    BOOST_CHECK(int(info.ntime()) == itsNTime);
    BOOST_CHECK(info.firstTime() == 100.0);
    BOOST_CHECK(info.timeInterval() == 5.0);
    BOOST_CHECK(int(info.nchanAvg()) == 1);
    BOOST_CHECK(int(info.ntimeAvg()) == 1);
    BOOST_CHECK(int(info.chanFreqs().size()) == itsNChan);
    BOOST_CHECK(int(info.chanWidths().size()) == itsNChan);
    BOOST_CHECK(info.msName().empty());
  }

  int itsCount;
  int itsNTime, itsNBl, itsNChan, itsNCorr;
};

// Test simple flagging with or without preflagged points.
void test1(int ntime, int nant, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("timewindow", "1");
  Step::ShPtr step2(new AOFlaggerStep(parset, ""));
  Step::ShPtr step3(new TestOutput(ntime, nant, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

// Test applyautocorr flagging with or without preflagged points.
void test2(int ntime, int nant, int nchan, int ncorr, bool flag) {
  // Create the steps.
  TestInput* in = new TestInput(ntime, nant, nchan, ncorr, flag);
  Step::ShPtr step1(in);
  ParameterSet parset;
  parset.add("timewindow", "4");
  parset.add("overlapmax", "1");
  Step::ShPtr step2(new AOFlaggerStep(parset, ""));
  Step::ShPtr step3(new TestOutput(ntime, nant, nchan, ncorr));
  dp3::steps::test::Execute({step1, step2, step3});
}

BOOST_AUTO_TEST_CASE(legacy_test1) {
  for (unsigned int i = 0; i < 2; ++i) {
    test1(10, 2, 32, 4, false);
    test1(10, 5, 32, 4, true);
  }
}

BOOST_AUTO_TEST_CASE(legacy_test2) {
  for (unsigned int i = 0; i < 2; ++i) {
    test2(4, 2, 8, 4, false);
    test2(10, 5, 32, 4, true);
    test2(8, 2, 8, 4, false);
    test2(14, 2, 8, 4, false);
  }
}

BOOST_AUTO_TEST_SUITE_END()
