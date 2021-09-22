// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "../../DDECal.h"
#include "../../InputStep.h"

BOOST_AUTO_TEST_SUITE(ddecal)

/// Helper forwarder to improve test failure output.
void TestShow(std::string expected,
              std::vector<std::pair<std::string, std::string>> parameters) {
  BOOST_CHECK_EQUAL(expected,
                    dp3::steps::test::Show<dp3::steps::DDECal>(
                        dp3::steps::test::CreateParameterSet(parameters)));
}

BOOST_AUTO_TEST_CASE(show_default) {
  TestShow(
      R"(DDECal prefix.
  mode (constraints):  diagonal
  algorithm:           directionsolve
  H5Parm:              tDDECal.MS/instrument.h5
  solint:              1
  nchan:               1
  directions:          [[center]]
  tolerance:           0.0001
  max iter:            50
  flag unconverged:    false
     diverged only:    false
  propagate solutions: false
       converged only: false
  detect stalling:     true
  step size:           0.2
  approximate fitter:  false
  only predict:        false
  subtract model:      false
Model steps for direction center
Predict
OnePredict prefix.
  sourcedb:           tDDECal.MS/sky
   number of patches: 1
   number of sources: 1
   all unpolarized:   true
   correct freq smearing: false
  apply beam:         false
  operation:          replace
  threads:            )" +
          std::to_string(aocommon::ThreadPool::NCPUs()) + R"(

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky"}});
}

BOOST_AUTO_TEST_CASE(show_modified) {
  TestShow(
      R"(DDECal prefix.
  mode (constraints):  diagonal
  algorithm:           hybrid
  H5Parm:              tDDECal.MS/instrument.h5
  solint:              42
  nchan:               44
  directions:          [[center]]
  min visib. ratio:    43.123
  tolerance:           1e-05
  max iter:            49
  flag unconverged:    true
     diverged only:    true
  propagate solutions: true
       converged only: true
  detect stalling:     true
  step size:           0.2
  coreconstraint:      45.123
  smoothnessconstraint:46.123
  approximate fitter:  false
  only predict:        true
  subtract model:      true
Model steps for direction center
Predict
OnePredict prefix.
  sourcedb:           tDDECal.MS/sky
   number of patches: 1
   number of sources: 1
   all unpolarized:   true
   correct freq smearing: false
  apply beam:         false
  operation:          replace
  threads:            )" +
          std::to_string(aocommon::ThreadPool::NCPUs()) + R"(

)",
      {{"msin", "tDDECal.MS"},
       {"prefix.directions", "[[center]]"},
       {"prefix.sourcedb", "tDDECal.MS/sky"},
       {"prefix.propagatesolutions", "true"},
       {"prefix.propagateconvergedonly", "true"},
       {"prefix.flagunconverged", "true"},
       {"prefix.flagdivergedonly", "true"},
       {"prefix.onlypredict", "true"},
       {"prefix.subtract", "true"},
       {"prefix.solveralgorithm", "hybrid"},
       {"prefix.solint", "42"},
       {"prefix.minvisratio", "43.123"},
       {"prefix.nchan", "44"},
       {"prefix.coreconstraint", "45.123"},
       {"prefix.smoothnessconstraint", "46.123"},
       {"prefix.maxiter", "49"}});
}

BOOST_AUTO_TEST_SUITE_END()
