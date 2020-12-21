// DDECal.cc: DPPP step class to do a direction dependent gain calibration
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// @author Tammo Jan Dijkema & André Offringa

#include "DDECal.h"

#include "../DPPP/DPBuffer.h"
#include "../DPPP/DPInfo.h"
#include "../DPPP/DPLogger.h"
#include "../DPPP/MSReader.h"
#include "../DPPP/Simulate.h"
#include "../DPPP/SourceDBUtil.h"
#include "../DPPP/Version.h"
#include "../DPPP/ColumnReader.h"
#include "../DPPP/DPRun.h"

#include "../IDGPredict/IDGPredict.h"

#include "DiagonalSolver.h"
#include "FullJonesSolver.h"
#include "RotationConstraint.h"
#include "RotationAndDiagonalConstraint.h"
#include "ScalarSolver.h"
#include "SmoothnessConstraint.h"
#include "TECConstraint.h"

#include <aocommon/matrix2x2.h>

#ifdef HAVE_ARMADILLO
#include "ScreenConstraint.h"
#endif

#include "../ParmDB/ParmDB.h"
#include "../ParmDB/ParmValue.h"
#include "../ParmDB/SourceDB.h"

#include "../Common/Memory.h"
#include "../Common/ParameterSet.h"
#include "../Common/StreamUtil.h"
#include "../Common/StringUtil.h"

#include <aocommon/threadpool.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/make_unique.hpp>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/OS/File.h>

#include <algorithm>
#include <cassert>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

using aocommon::FitsReader;
using namespace casacore;
using namespace DP3::BBS;
using schaapcommon::h5parm::SolTab;

namespace DP3 {
namespace DPPP {

DDECal::DDECal(DPInput* input, const ParameterSet& parset, const string& prefix)
    : itsInput(input),
      itsName(prefix),
      itsAvgTime(0),
      itsSols(),
      itsH5ParmName(parset.getString(
          prefix + "h5parm", parset.getString("msin") + "/instrument.h5")),
      itsH5Parm(itsH5ParmName, true),
      itsPropagateSolutions(
          parset.getBool(prefix + "propagatesolutions", false)),
      itsPropagateConvergedOnly(
          parset.getBool(prefix + "propagateconvergedonly", false)),
      itsFlagUnconverged(parset.getBool(prefix + "flagunconverged", false)),
      itsFlagDivergedOnly(parset.getBool(prefix + "flagdivergedonly", false)),
      itsOnlyPredict(parset.getBool(prefix + "onlypredict", false)),
      itsTimeStep(0),
      itsSolInt(parset.getInt(prefix + "solint", 1)),
      itsSolIntCount(1),
      itsNSolInts(0),
      itsMinVisRatio(parset.getDouble(prefix + "minvisratio", 0.0)),
      itsStepInSolInt(0),
      itsBufferedSolInts(0),
      itsNChan(parset.getInt(prefix + "nchan", 1)),
      itsUVWFlagStep(input, parset, prefix),
      itsCoreConstraint(parset.getDouble(prefix + "coreconstraint", 0.0)),
      itsAntennaConstraint(),
      itsSmoothnessConstraint(
          parset.getDouble(prefix + "smoothnessconstraint", 0.0)),
      itsScreenCoreConstraint(
          parset.getDouble(prefix + "tecscreen.coreconstraint", 0.0)),
      itsPolsInSolutions(1),
      itsApproximateTEC(false),
      itsSubtract(parset.getBool(prefix + "subtract", false)),
      itsStatFilename(parset.getString(prefix + "statfilename", "")) {
  stringstream ss;
  ss << parset;
  itsParsetString = ss.str();

  if (!itsStatFilename.empty())
    itsStatStream = boost::make_unique<std::ofstream>(itsStatFilename);

  // Read the antennaconstraint list
  std::vector<std::string> antConstraintList = parset.getStringVector(
      prefix + "antennaconstraint", std::vector<std::string>());
  if (!antConstraintList.empty()) {
    for (const std::string& antSetStr : antConstraintList) {
      ParameterValue antSetParam(antSetStr);
      std::vector<std::string> list = antSetParam.getStringVector();
      itsAntennaConstraint.emplace_back(list.begin(), list.end());
      // By doing this check after inserting in the set, duplicate antenna names
      // will be removed.
      if (itsAntennaConstraint.back().size() == 1)
        throw std::runtime_error(
            "Error: antennaconstraint given that should constrain a group of "
            "antennas with one antenna in it. This does not make sense (did "
            "you forget to use two square brackets? [[ ant1, ant2 ]] )");
    }
  }

  itsMode = GainCal::stringToCalType(
      boost::to_lower_copy(parset.getString(prefix + "mode", "complexgain")));

  initializeSolver(parset, prefix);

  itsSolver->SetMaxIterations(parset.getInt(prefix + "maxiter", 50));
  double tolerance = parset.getDouble(prefix + "tolerance", 1.e-4);
  itsSolver->SetAccuracy(tolerance);
  itsSolver->SetConstraintAccuracy(
      parset.getDouble(prefix + "approxtolerance", tolerance * 10.0));
  itsSolver->SetStepSize(parset.getDouble(prefix + "stepsize", 0.2));
  itsSolver->SetDetectStalling(parset.getBool(prefix + "detectstalling", true));

  // Initialize steps
  initializeColumnReaders(parset, prefix);
  initializeIDG(parset, prefix);
  initializePredictSteps(parset, prefix);

  if (itsDirections.size() == 0)
    throw std::runtime_error(
        "DDECal initialized with 0 directions: something is wrong with your "
        "parset or your sourcedb");
}

DDECal::~DDECal() {}

DPStep::ShPtr DDECal::makeStep(DPInput* input, const ParameterSet& parset,
                               const std::string& prefix) {
  return std::make_shared<DDECal>(input, parset, prefix);
}

void DDECal::initializeSolver(const ParameterSet& parset,
                              const string& prefix) {
  if (itsCoreConstraint != 0.0 || !itsAntennaConstraint.empty()) {
    itsConstraints.push_back(boost::make_unique<AntennaConstraint>());
  }
  if (itsSmoothnessConstraint != 0.0) {
    itsConstraints.emplace_back(
        new SmoothnessConstraint(itsSmoothnessConstraint));
  }
  switch (itsMode) {
    case GainCal::SCALARCOMPLEXGAIN:
      itsSolver = boost::make_unique<ScalarSolver>();
      // no constraints
      itsSolver->SetPhaseOnly(false);
      itsPolsInSolutions = 1;
      break;
    case GainCal::SCALARPHASE:
      itsSolver = boost::make_unique<ScalarSolver>();
      itsConstraints.push_back(boost::make_unique<PhaseOnlyConstraint>());
      itsSolver->SetPhaseOnly(true);
      itsPolsInSolutions = 1;
      break;
    case GainCal::SCALARAMPLITUDE:
      itsSolver = boost::make_unique<ScalarSolver>();
      itsConstraints.push_back(boost::make_unique<AmplitudeOnlyConstraint>());
      itsSolver->SetPhaseOnly(false);
      itsPolsInSolutions = 1;
      break;
    case GainCal::DIAGONAL:
      itsSolver = boost::make_unique<DiagonalSolver>();
      // no constraints
      itsSolver->SetPhaseOnly(false);
      itsPolsInSolutions = 2;
      break;
    case GainCal::DIAGONALPHASE:
      itsSolver = boost::make_unique<DiagonalSolver>();
      itsConstraints.push_back(boost::make_unique<PhaseOnlyConstraint>());
      itsSolver->SetPhaseOnly(true);
      itsPolsInSolutions = 2;
      break;
    case GainCal::DIAGONALAMPLITUDE:
      itsSolver = boost::make_unique<DiagonalSolver>();
      itsConstraints.push_back(boost::make_unique<AmplitudeOnlyConstraint>());
      itsSolver->SetPhaseOnly(false);
      itsPolsInSolutions = 2;
      break;
    case GainCal::FULLJONES:
      itsSolver = boost::make_unique<FullJonesSolver>();
      // no constraints
      itsSolver->SetPhaseOnly(false);
      itsPolsInSolutions = 4;
      break;
    case GainCal::TEC:
    case GainCal::TECANDPHASE: {
      itsSolver = boost::make_unique<ScalarSolver>();
      const auto tecMode = (itsMode == GainCal::TEC)
                               ? TECConstraint::TECOnlyMode
                               : TECConstraint::TECAndCommonScalarMode;
      std::unique_ptr<TECConstraint> constraintPtr;

      itsApproximateTEC = parset.getBool(prefix + "approximatetec", false);
      if (itsApproximateTEC) {
        const int iters = parset.getInt(prefix + "maxapproxiter",
                                        itsSolver->GetMaxIterations() / 2);
        const int chunksize = parset.getInt(prefix + "approxchunksize", 0);
        auto approxConstraint =
            boost::make_unique<ApproximateTECConstraint>(tecMode);
        approxConstraint->SetMaxApproximatingIterations(iters);
        approxConstraint->SetFittingChunkSize(chunksize);
        constraintPtr = std::move(approxConstraint);
      } else {
        constraintPtr = boost::make_unique<TECConstraint>(tecMode);
      }
      constraintPtr->setDoPhaseReference(
          parset.getBool(prefix + "phasereference", true));
      itsConstraints.push_back(std::move(constraintPtr));
      itsSolver->SetPhaseOnly(true);
      itsPolsInSolutions = 1;
      break;
    }
    case GainCal::TECSCREEN:
#ifdef HAVE_ARMADILLO
      itsSolver = boost::make_unique<ScalarSolver>();
      itsConstraints.push_back(
          boost::make_unique<ScreenConstraint>(parset, prefix + "tecscreen."));
      itsSolver->SetPhaseOnly(true);
      itsPolsInSolutions = 1;
#else
      throw std::runtime_error(
          "Can not use TEC screen: Armadillo is not available. Recompile DP3 "
          "with Armadillo.");
#endif
      break;
    case GainCal::ROTATIONANDDIAGONAL: {
      itsSolver = boost::make_unique<FullJonesSolver>();
      auto constraintPtr = boost::make_unique<RotationAndDiagonalConstraint>();
      constraintPtr->SetDoRotationReference(
          parset.getBool(prefix + "rotationreference", false));
      itsConstraints.push_back(std::move(constraintPtr));
      itsPolsInSolutions = 4;
      break;
    }
    case GainCal::ROTATION:
      itsSolver = boost::make_unique<FullJonesSolver>();
      itsConstraints.push_back(boost::make_unique<RotationConstraint>());
      itsPolsInSolutions = 4;
      break;
    default:
      throw std::runtime_error("Unexpected mode: " +
                               GainCal::calTypeToString(itsMode));
  }
}

void DDECal::initializeColumnReaders(const ParameterSet& parset,
                                     const string& prefix) {
  std::vector<std::string> cols = parset.getStringVector(
      prefix + "modeldatacolumns", std::vector<std::string>());

  // The statement below allows DDECal to be backwards compatible, e.g.
  // DPPP msin.modelcolumn=MY_MODEL_DATA ddecal.usemodelcolumn=true
  // msin=tDDECal.MS msout=.
  if (cols.size() == 0 && parset.getBool(prefix + "usemodelcolumn", false)) {
    cols.push_back(parset.getString("msin.modelcolumn", "MODEL_DATA"));
  }
  for (string& col : cols) {
    if (cols.size() == 1) {
      itsDirections.emplace_back(1, "pointing");
    } else {
      itsDirections.emplace_back(1, col);
    }
    itsSteps.push_back(
        std::make_shared<ColumnReader>(*itsInput, parset, prefix, col));
    setModelNextSteps(itsSteps.back(), col, parset, prefix);
  }
}

void DDECal::initializeIDG(const ParameterSet& parset, const string& prefix) {
  // TODO it would be nicer to get a new method in the DS9 reader to first get
  // names of directions and pass that to idgpredict. It will then read it
  // itself instead of DDECal having to do everything. It is better to do it all
  // in IDGPredict, so we can also make it
  std::string regionFilename = parset.getString(prefix + "idg.regions", "");
  std::vector<std::string> imageFilenames =
      parset.getStringVector(prefix + "idg.images", std::vector<string>());

  if (regionFilename.empty() && imageFilenames.empty()) {
    return;
  }

  std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
      readers = IDGPredict::GetReaders(imageFilenames);
  std::vector<Facet> facets =
      IDGPredict::GetFacets(regionFilename, readers.first.front());

  for (size_t i = 0; i < facets.size(); ++i) {
    std::string dir_name = "dir" + std::to_string(i);
    if (!facets[i].Direction().empty()) {
      dir_name = facets[i].Direction();
    }
    itsDirections.emplace_back(1, dir_name);

    itsSteps.push_back(std::make_shared<IDGPredict>(
        *itsInput, parset, prefix, readers, std::vector<Facet>{facets[i]}));
    setModelNextSteps(itsSteps.back(), facets[i].Direction(), parset, prefix);
  }
}

void DDECal::initializePredictSteps(const ParameterSet& parset,
                                    const string& prefix) {
  std::vector<string> strDirections =
      parset.getStringVector(prefix + "directions", std::vector<string>());
  size_t start = itsDirections.size();
  string sourceDBName = parset.getString(prefix + "sourcedb", "");
  // Default directions are all patches
  if (strDirections.empty() && !sourceDBName.empty()) {
    BBS::SourceDB sourceDB(BBS::ParmDBMeta("", sourceDBName), false);
    std::vector<string> patchNames =
        makePatchList(sourceDB, std::vector<string>());
    for (const string& patch : patchNames) {
      itsDirections.emplace_back(1, patch);
    }
  } else {
    for (const string& direction : strDirections) {
      ParameterValue dirStr(direction);
      itsDirections.emplace_back(dirStr.getStringVector());
    }
  }

  for (size_t dir = start; dir < itsDirections.size(); ++dir) {
    itsSteps.push_back(std::make_shared<Predict>(itsInput, parset, prefix,
                                                 itsDirections[dir]));
    setModelNextSteps(itsSteps.back(), itsDirections[dir][0], parset, prefix);
  }
}

void DDECal::setModelNextSteps(std::shared_ptr<DPStep> step,
                               const std::string direction,
                               const ParameterSet& parset,
                               const string prefix) {
  std::string current_steps = parset.getString("steps");
  std::string model_next_steps =
      parset.getString(prefix + "modelnextsteps", "[]");
  if (parset.isDefined(prefix + "modelnextsteps." + direction)) {
    model_next_steps = parset.getString(prefix + "modelnextsteps." + direction);
  }

  // Make a shallow copy to work around constness of parset
  ParameterSet parset_new(parset);
  parset_new.replace("steps", model_next_steps);

  auto first_step = DPRun::makeSteps(parset_new, "", itsInput, false);

  if (first_step) {
    step->setNextStep(first_step);
  }

  // Revert the changes made to the shallow copy
  parset_new.replace("steps", current_steps);
}

void DDECal::updateInfo(const DPInfo& infoIn) {
  DPStep::updateInfo(infoIn);
  info().setNeedVisData();
  if (itsSubtract) info().setWriteData();

  const size_t nDir = itsDirections.size();

  itsUVWFlagStep.updateInfo(infoIn);
  itsThreadPool =
      boost::make_unique<aocommon::ThreadPool>(getInfo().nThreads());

  // Update info for substeps and set other required parameters
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    itsSteps[dir]->setInfo(infoIn);

    if (auto s = std::dynamic_pointer_cast<Predict>(itsSteps[dir])) {
      s->setThreadData(*itsThreadPool, itsMeasuresMutex);
    } else if (auto s = std::dynamic_pointer_cast<IDGPredict>(itsSteps[dir])) {
      itsSolIntCount = std::max(
          itsSolIntCount, s->GetBufferSize() / itsSteps.size() / itsSolInt);
      // We increment by one so the IDGPredict will not flush in its process
      s->SetBufferSize(itsSolInt * itsSolIntCount + 1);
    } else if (!std::dynamic_pointer_cast<ColumnReader>(itsSteps[dir])) {
      throw std::runtime_error("DDECal received an invalid first model step");
    }
  }

  itsSolver->SetNThreads(getInfo().nThreads());

  if (itsSolInt == 0) {
    itsSolInt = info().ntime();
  }

  for (std::unique_ptr<Constraint>& constraint : itsConstraints) {
    constraint->SetNThreads(getInfo().nThreads());
    itsSolver->AddConstraint(*constraint);
  }

  itsDataResultStep = std::make_shared<ResultStep>();
  itsUVWFlagStep.setNextStep(itsDataResultStep);

  // Each step should be coupled to a resultstep
  itsResultSteps.resize(itsSteps.size());
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    itsResultSteps[dir] =
        std::make_shared<MultiResultStep>(itsSolInt * itsSolIntCount);

    // Add the resultstep to the end of the model next steps
    std::shared_ptr<DPStep> step = itsSteps[dir];
    while (step->getNextStep()) {
      step = step->getNextStep();
    }
    step->setNextStep(itsResultSteps[dir]);
  }

  if (itsNChan == 0 || itsNChan > info().nchan()) {
    itsNChan = info().nchan();
  }

  // Convert from casacore::Vector to std::vector, pass only used antennas to
  // multidirsolver
  std::vector<int> ant1(info().getAnt1().size());
  std::vector<int> ant2(info().getAnt2().size());
  for (unsigned int i = 0; i < ant1.size(); ++i) {
    ant1[i] = info().antennaMap()[info().getAnt1()[i]];
    ant2[i] = info().antennaMap()[info().getAnt2()[i]];
  }

  // Fill antenna info in H5Parm, need to convert from casa types to std types
  // Fill in metadata for all antennas, also those that may be filtered out.
  std::vector<std::string> antennaNames(info().antennaNames().size());
  std::vector<std::vector<double>> antennaPos(info().antennaPos().size());
  for (unsigned int i = 0; i < info().antennaNames().size(); ++i) {
    antennaNames[i] = info().antennaNames()[i];
    casacore::Quantum<casacore::Vector<double>> pos =
        info().antennaPos()[i].get("m");
    antennaPos[i].resize(3);
    antennaPos[i][0] = pos.getValue()[0];
    antennaPos[i][1] = pos.getValue()[1];
    antennaPos[i][2] = pos.getValue()[2];
  }

  itsH5Parm.AddAntennas(antennaNames, antennaPos);

  std::vector<std::pair<double, double>> sourcePositions(itsDirections.size());

  for (unsigned int i = 0; i < itsSteps.size(); ++i) {
    if (auto s = std::dynamic_pointer_cast<Predict>(itsSteps[i])) {
      sourcePositions[i] = s->getFirstDirection();
    } else if (auto s = std::dynamic_pointer_cast<IDGPredict>(itsSteps[i])) {
      // We can take the front element of an IDG step since it only contains 1.
      sourcePositions[i] = s->GetFirstDirection();
    } else if (auto s = std::dynamic_pointer_cast<ColumnReader>(itsSteps[i])) {
      MDirection dirJ2000(
          MDirection::Convert(infoIn.phaseCenter(), MDirection::J2000)());
      Quantum<Vector<Double>> angles = dirJ2000.getAngle();
      sourcePositions[i] = std::pair<double, double>(angles.getBaseValue()[0],
                                                     angles.getBaseValue()[1]);
    }
  }

  itsH5Parm.AddSources(getDirectionNames(), sourcePositions);

  size_t nSolTimes = (info().ntime() + itsSolInt - 1) / itsSolInt;
  size_t nChannelBlocks = info().nchan() / itsNChan;
  itsSols.resize(nSolTimes);
  itsNIter.resize(nSolTimes);
  itsNApproxIter.resize(nSolTimes);
  itsConstraintSols.resize(nSolTimes);
  itsVisInInterval.assign(nChannelBlocks, std::pair<size_t, size_t>(0, 0));

  itsChanBlockStart.resize(nChannelBlocks + 1);
  itsChanBlockFreqs.resize(nChannelBlocks);
  for (size_t chBlock = 0; chBlock != nChannelBlocks; ++chBlock) {
    const size_t channelIndexStart = chBlock * info().nchan() / nChannelBlocks,
                 channelIndexEnd =
                     (chBlock + 1) * info().nchan() / nChannelBlocks,
                 curChannelBlockSize = channelIndexEnd - channelIndexStart;
    double meanfreq =
        std::accumulate(info().chanFreqs().data() + channelIndexStart,
                        info().chanFreqs().data() + channelIndexEnd, 0.0) /
        curChannelBlockSize;
    itsChanBlockStart[chBlock] = channelIndexStart;
    itsChanBlockFreqs[chBlock] = meanfreq;
  }
  itsChanBlockStart.back() = info().nchan();

  itsWeightsPerAntenna.assign(
      itsChanBlockFreqs.size() * info().antennaUsed().size(), 0.0);

  for (unsigned int i = 0; i < itsConstraints.size(); ++i) {
    // Initialize the constraint with some common metadata
    itsConstraints[i]->InitializeDimensions(
        info().antennaNames().size(), itsDirections.size(), nChannelBlocks);

    // Different constraints need different information. Determine if the
    // constraint is of a type that needs more information, and if so initialize
    // the constraint.
    AntennaConstraint* antConstraint =
        dynamic_cast<AntennaConstraint*>(itsConstraints[i].get());
    if (antConstraint != nullptr) {
      if (itsAntennaConstraint.empty()) {
        // Set the antenna constraint to all stations within certain distance
        // specified by 'coreconstraint' parameter.
        // Take antenna with index 0 as reference station
        double refX = antennaPos[0][0], refY = antennaPos[0][1],
               refZ = antennaPos[0][2];
        std::vector<std::set<size_t>> antConstraintList(1);
        std::set<size_t>& coreAntennaIndices = antConstraintList.front();
        const double coreDistSq = itsCoreConstraint * itsCoreConstraint;
        for (size_t ant = 0; ant != antennaPos.size(); ++ant) {
          double dx = refX - antennaPos[ant][0], dy = refY - antennaPos[ant][1],
                 dz = refZ - antennaPos[ant][2],
                 distSq = dx * dx + dy * dy + dz * dz;
          if (distSq <= coreDistSq) coreAntennaIndices.insert(ant);
        }
        antConstraint->initialize(std::move(antConstraintList));
      } else {
        // Set the antenna constraint to a list of stations indices that
        // are to be kept the same during the solve.
        const casacore::Vector<casacore::String>& antNames =
            info().antennaNames();
        std::vector<std::string> antNamesStl(
            antNames.begin(),
            antNames.end());  // casacore vector doesn't support find properly
        std::vector<std::set<size_t>> constraintList;
        for (const std::set<std::string>& constraintNameSet :
             itsAntennaConstraint) {
          constraintList.emplace_back();
          for (const std::string& constraintName : constraintNameSet) {
            auto iter = std::find(antNamesStl.begin(), antNamesStl.end(),
                                  constraintName);
            if (iter != antNamesStl.end())
              constraintList.back().insert(
                  info().antennaMap()[iter - antNamesStl.begin()]);
          }
          if (constraintList.back().size() <= 1)
            throw std::runtime_error(
                "Error in antenna constraint: at least two antennas expected");
        }
        antConstraint->initialize(std::move(constraintList));
      }
    }

#ifdef HAVE_ARMADILLO
    ScreenConstraint* screenConstraint =
        dynamic_cast<ScreenConstraint*>(itsConstraints[i].get());
    if (screenConstraint != 0) {
      screenConstraint->initialize(&(itsChanBlockFreqs[0]));
      screenConstraint->setAntennaPositions(antennaPos);
      screenConstraint->setDirections(sourcePositions);
      screenConstraint->initPiercePoints();
      double refX = antennaPos[i][0], refY = antennaPos[i][1],
             refZ = antennaPos[i][2];
      std::vector<size_t> coreAntennaIndices;
      std::vector<size_t> otherAntennaIndices;
      const double coreDistSq =
          itsScreenCoreConstraint * itsScreenCoreConstraint;
      for (size_t ant = 0; ant != antennaPos.size(); ++ant) {
        double dx = refX - antennaPos[ant][0], dy = refY - antennaPos[ant][1],
               dz = refZ - antennaPos[ant][2],
               distSq = dx * dx + dy * dy + dz * dz;
        if (distSq <= coreDistSq)
          coreAntennaIndices.emplace_back(info().antennaMap()[ant]);
        else
          otherAntennaIndices.emplace_back(info().antennaMap()[ant]);
      }
      screenConstraint->setCoreAntennas(coreAntennaIndices);
      screenConstraint->setOtherAntennas(otherAntennaIndices);
    }
#endif

    TECConstraintBase* tecConstraint =
        dynamic_cast<TECConstraintBase*>(itsConstraints[i].get());
    if (tecConstraint != nullptr) {
      tecConstraint->initialize(&itsChanBlockFreqs[0]);
    }
    SmoothnessConstraint* sConstraint =
        dynamic_cast<SmoothnessConstraint*>(itsConstraints[i].get());
    if (sConstraint != nullptr) {
      sConstraint->Initialize(&itsChanBlockFreqs[0]);
    }
  }

  unsigned int nSt = info().antennaUsed().size();
  // Give renumbered antennas to multidirsolver

  itsSolver->Initialize(nSt, nDir, info().nchan(), nChannelBlocks, ant1, ant2);

  for (unsigned int i = 0; i < nSolTimes; ++i) {
    itsSols[i].resize(nChannelBlocks);
  }
}

void DDECal::show(std::ostream& os) const {
  os << "DDECal " << itsName << '\n'
     << "  H5Parm:              " << itsH5ParmName << '\n'
     << "  solint:              " << itsSolInt << '\n'
     << "  nchan:               " << itsNChan << '\n'
     << "  directions:          " << itsDirections << '\n';
  if (itsMinVisRatio != 0.0) {
    os << "  min visib. ratio:    " << itsMinVisRatio << '\n';
  }
  os << "  tolerance:           " << itsSolver->GetAccuracy() << '\n'
     << "  max iter:            " << itsSolver->GetMaxIterations() << '\n'
     << "  flag unconverged:    " << std::boolalpha << itsFlagUnconverged
     << '\n'
     << "     diverged only:    " << std::boolalpha << itsFlagDivergedOnly
     << '\n'
     << "  propagate solutions: " << std::boolalpha << itsPropagateSolutions
     << '\n'
     << "       converged only: " << std::boolalpha << itsPropagateConvergedOnly
     << '\n'
     << "  detect stalling:     " << std::boolalpha
     << itsSolver->GetDetectStalling() << '\n'
     << "  step size:           " << itsSolver->GetStepSize() << '\n'
     << "  mode (constraints):  " << GainCal::calTypeToString(itsMode) << '\n';
  if (!itsAntennaConstraint.empty())
    os << "  antennaconstraint:   " << itsAntennaConstraint << '\n';
  if (itsCoreConstraint != 0.0)
    os << "  coreconstraint:      " << itsCoreConstraint << '\n';
  if (itsSmoothnessConstraint != 0.0)
    os << "  smoothnessconstraint:" << itsSmoothnessConstraint << '\n';
  os << "  approximate fitter:  " << itsApproximateTEC << '\n'
     << "  only predict:        " << itsOnlyPredict << '\n'
     << "  subtract model:      " << itsSubtract << '\n';
  for (unsigned int i = 0; i < itsSteps.size(); ++i) {
    std::shared_ptr<DPStep> step = itsSteps[i];
    os << "Model steps for direction " << itsDirections[i][0] << '\n';
    do {
      step->show(os);
    } while (step = step->getNextStep());
    os << '\n';
  }
  itsUVWFlagStep.show(os);
}

void DDECal::showTimings(std::ostream& os, double duration) const {
  double totaltime = itsTimer.getElapsed();
  os << "  ";
  FlagCounter::showPerc1(os, itsTimer.getElapsed(), duration);
  os << " DDECal " << itsName << endl;

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerPredict.getElapsed(), totaltime);
  os << " of it spent in predict" << endl;

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerSolve.getElapsed(), totaltime);
  os << " of it spent in estimating gains and computing residuals" << endl;

  itsSolver->GetTimings(os, itsTimerSolve.getElapsed());

  os << "          ";
  FlagCounter::showPerc1(os, itsTimerWrite.getElapsed(), totaltime);
  os << " of it spent in writing gain solutions to disk" << endl;

  os << "          ";
  os << "Substeps taken:" << std::endl;
  for (auto& step : itsSteps) {
    os << "          ";
    step->showTimings(os, duration);
  }

  os << "Iterations taken: [";
  for (unsigned int i = 0; i < itsNIter.size() - 1; ++i) {
    os << itsNIter[i];
    if (itsNApproxIter[i] != 0) os << '|' << itsNApproxIter[i];
    os << ",";
  }
  os << itsNIter[itsNIter.size() - 1];
  if (itsNApproxIter[itsNIter.size() - 1] != 0)
    os << '|' << itsNApproxIter[itsNIter.size() - 1];
  os << "]" << endl;
}

void DDECal::initializeScalarSolutions(size_t bufferIndex) {
  if (sol_ints_[bufferIndex].NSolution() > 0 && itsPropagateSolutions) {
    if (itsNIter[sol_ints_[bufferIndex].NSolution() - 1] >
            itsSolver->GetMaxIterations() &&
        itsPropagateConvergedOnly) {
      // initialize solutions with 1.
      size_t n = itsDirections.size() * info().antennaUsed().size();
      for (std::vector<DComplex>& solvec :
           itsSols[sol_ints_[bufferIndex].NSolution()]) {
        solvec.assign(n, 1.0);
      }
    } else {
      // initialize solutions with those of the previous step
      itsSols[sol_ints_[bufferIndex].NSolution()] =
          itsSols[sol_ints_[bufferIndex].NSolution() - 1];
    }
  } else {
    // initialize solutions with 1.
    size_t n = itsDirections.size() * info().antennaUsed().size();
    for (std::vector<DComplex>& solvec :
         itsSols[sol_ints_[bufferIndex].NSolution()]) {
      solvec.assign(n, 1.0);
    }
  }
}

void DDECal::initializeFullMatrixSolutions(size_t bufferIndex) {
  if (sol_ints_[bufferIndex].NSolution() > 0 && itsPropagateSolutions) {
    if (itsNIter[sol_ints_[bufferIndex].NSolution() - 1] >
            itsSolver->GetMaxIterations() &&
        itsPropagateConvergedOnly) {
      // initialize solutions with unity matrix [1 0 ; 0 1].
      size_t n = itsDirections.size() * info().antennaUsed().size();
      for (std::vector<DComplex>& solvec :
           itsSols[sol_ints_[bufferIndex].NSolution()]) {
        solvec.resize(n * 4);
        for (size_t i = 0; i != n; ++i) {
          solvec[i * 4 + 0] = 1.0;
          solvec[i * 4 + 1] = 0.0;
          solvec[i * 4 + 2] = 0.0;
          solvec[i * 4 + 3] = 1.0;
        }
      }
    } else {
      // initialize solutions with those of the previous step
      itsSols[sol_ints_[bufferIndex].NSolution()] =
          itsSols[sol_ints_[bufferIndex].NSolution() - 1];
    }
  } else {
    // initialize solutions with unity matrix [1 0 ; 0 1].
    size_t n = itsDirections.size() * info().antennaUsed().size();
    for (std::vector<DComplex>& solvec :
         itsSols[sol_ints_[bufferIndex].NSolution()]) {
      solvec.resize(n * 4);
      for (size_t i = 0; i != n; ++i) {
        solvec[i * 4 + 0] = 1.0;
        solvec[i * 4 + 1] = 0.0;
        solvec[i * 4 + 2] = 0.0;
        solvec[i * 4 + 3] = 1.0;
      }
    }
  }
}

std::vector<std::string> DDECal::getDirectionNames() {
  std::vector<std::string> res;

  for (const std::vector<std::string>& dir : itsDirections) {
    std::stringstream ss;
    ss << dir;
    res.emplace_back(ss.str());
  }

  return res;
}

void DDECal::flagChannelBlock(size_t cbIndex, size_t bufferIndex) {
  const size_t nBl = info().nbaselines(), nCh = info().nchan(),
               nChanBlocks = itsChanBlockFreqs.size();
  // Set the antenna-based weights to zero
  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t ant1 = info().antennaMap()[info().getAnt1()[bl]];
    size_t ant2 = info().antennaMap()[info().getAnt2()[bl]];
    for (size_t ch = itsChanBlockStart[cbIndex];
         ch != itsChanBlockStart[cbIndex + 1]; ++ch) {
      itsWeightsPerAntenna[ant1 * nChanBlocks + cbIndex] = 0.0;
      itsWeightsPerAntenna[ant2 * nChanBlocks + cbIndex] = 0.0;
    }
  }
  // Set the visibility weights to zero
  for (size_t step = 0; step != sol_ints_[bufferIndex].Size(); ++step) {
    for (size_t bl = 0; bl < nBl; ++bl) {
      for (size_t chcr = itsChanBlockStart[cbIndex] * 4;
           chcr != itsChanBlockStart[cbIndex + 1] * 4; ++chcr) {
        const size_t index = bl * nCh * 4 + chcr;
        sol_ints_[bufferIndex].WeightPtrs()[step][index] = 0;
      }
    }
  }
}

void DDECal::checkMinimumVisibilities(size_t bufferIndex) {
  for (size_t cb = 0; cb != itsChanBlockFreqs.size(); ++cb) {
    double fraction =
        double(itsVisInInterval[cb].first) / itsVisInInterval[cb].second;
    if (fraction < itsMinVisRatio) flagChannelBlock(cb, bufferIndex);
  }
}

void DDECal::doSolve() {
  for (size_t dir = 0; dir < itsSteps.size(); ++dir) {
    if (auto s = dynamic_cast<IDGPredict*>(itsSteps[dir].get())) {
      itsTimerPredict.start();
      s->flush();
      itsTimerPredict.stop();
    }
    for (size_t i = 0; i < itsResultSteps[dir].get()->size(); ++i) {
      const size_t sol_int = i / itsSolInt;
      const size_t timestep = i % itsSolInt;

      sol_ints_[sol_int].ModelDataPtrs()[timestep][dir] =
          itsResultSteps[dir]->get()[i].getData().data();
    }
  }

  for (size_t i = 0; i < sol_ints_.size(); ++i) {
    // When the model data is subtracted after calibration, the model data
    // needs to be stored before solving, because the solver modifies it.
    // This is done conditionally to prevent using memory when it is
    // not required (the model data can be large).
    if (itsSubtract || itsOnlyPredict) {
      storeModelData(sol_ints_[i].ModelDataPtrs());
    }

    SolverBase::SolveResult solveResult;
    if (!itsOnlyPredict) {
      checkMinimumVisibilities(i);

      for (std::unique_ptr<Constraint>& constraint : itsConstraints) {
        constraint->SetWeights(itsWeightsPerAntenna);
      }

      if (itsPolsInSolutions == 2 || itsPolsInSolutions == 4)
        initializeFullMatrixSolutions(i);
      else
        initializeScalarSolutions(i);

      itsTimerSolve.start();

      // TODO to be done polymorphically once the solvers have been refactored
      if (dynamic_cast<DiagonalSolver*>(itsSolver.get())) {
        // Temporary fix: convert solutions from full Jones matrices to diagonal
        std::vector<std::vector<casacore::DComplex>>& full_solutions =
            itsSols[sol_ints_[i].NSolution()];
        std::vector<std::vector<std::complex<double>>> diagonals(
            full_solutions.size());
        for (size_t ch_block = 0; ch_block != diagonals.size(); ++ch_block) {
          diagonals[ch_block].reserve(full_solutions[ch_block].size() / 2);
          for (size_t s = 0; s != full_solutions[ch_block].size() / 4; ++s) {
            diagonals[ch_block].push_back(
                itsSols[sol_ints_[i].NSolution()][ch_block][s * 4]);
            diagonals[ch_block].push_back(
                itsSols[sol_ints_[i].NSolution()][ch_block][s * 4 + 3]);
          }
        }

        solveResult =
            itsSolver->Solve(sol_ints_[i].DataPtrs(), sol_ints_[i].WeightPtrs(),
                             std::move(sol_ints_[i].ModelDataPtrs()), diagonals,
                             itsAvgTime / itsSolInt, itsStatStream.get());

        // Temporary fix: extend solutions from diagonal to full Jones matrices
        for (size_t chBlock = 0; chBlock != full_solutions.size(); ++chBlock) {
          for (size_t s = 0; s != full_solutions[chBlock].size() / 4; ++s) {
            full_solutions[chBlock][s * 4] = diagonals[chBlock][s * 2];
            full_solutions[chBlock][s * 4 + 1] = 0.0;
            full_solutions[chBlock][s * 4 + 2] = 0.0;
            full_solutions[chBlock][s * 4 + 3] = diagonals[chBlock][s * 2 + 1];
          }
        }
      } else {
        solveResult =
            itsSolver->Solve(sol_ints_[i].DataPtrs(), sol_ints_[i].WeightPtrs(),
                             std::move(sol_ints_[i].ModelDataPtrs()),
                             itsSols[sol_ints_[i].NSolution()],
                             itsAvgTime / itsSolInt, itsStatStream.get());
      }
      itsTimerSolve.stop();

      itsNIter[sol_ints_[i].NSolution()] = solveResult.iterations;
      itsNApproxIter[sol_ints_[i].NSolution()] =
          solveResult.constraint_iterations;
    }

    if (itsSubtract || itsOnlyPredict) {
      subtractCorrectedModel(itsPolsInSolutions != 1, i);
    }

    // Check for nonconvergence and flag if desired. Unconverged solutions are
    // identified by the number of iterations being one more than the max
    // allowed number
    if (solveResult.iterations > itsSolver->GetMaxIterations() &&
        itsFlagUnconverged) {
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          if (itsFlagDivergedOnly) {
            // Set weights with negative values (indicating unconverged
            // solutions that diverged) to zero (all other unconverged
            // solutions remain unflagged)
            for (double& weight : result.weights) {
              if (weight < 0.) weight = 0.;
            }
          } else {
            // Set all weights to zero
            result.weights.assign(result.weights.size(), 0.);
          }
        }
      }
    } else {
      // Set any negative weights (indicating unconverged solutions that
      // diverged) to one (all other unconverged solutions are unflagged
      // already)
      for (auto& constraint_results : solveResult.results) {
        for (auto& result : constraint_results) {
          for (double& weight : result.weights) {
            if (weight < 0.0) weight = 1.0;
          }
        }
      }
    }

    // Store constraint solutions if any constaint has a non-empty result
    bool someConstraintHasResult = false;
    for (const auto& constraint_results : solveResult.results) {
      if (!constraint_results.empty()) {
        someConstraintHasResult = true;
        break;
      }
    }
    if (someConstraintHasResult) {
      itsConstraintSols[sol_ints_[i].NSolution()] = solveResult.results;
    }
  }

  itsTimer.stop();

  for (size_t i = 0; i < sol_ints_.size(); ++i) {
    sol_ints_[i].RestoreFlagsAndWeights();
    for (size_t step = 0; step < sol_ints_[i].Size(); ++step) {
      // Push data (possibly changed) to next step
      getNextStep()->process(sol_ints_[i][step]);
    }
  }

  itsTimer.start();
}

bool DDECal::process(const DPBuffer& bufin) {
  itsTimer.start();

  // Create a new solution interval if needed
  if (itsStepInSolInt == 0) {
    sol_ints_.emplace_back(itsInput, itsNSolInts, itsSolInt,
                           itsDirections.size(), itsTimer);
  }

  sol_ints_[itsBufferedSolInts].CopyBuffer(bufin);
  doPrepare(sol_ints_.back()[itsStepInSolInt], itsBufferedSolInts,
            itsStepInSolInt);

  ++itsStepInSolInt;

  if (itsStepInSolInt == itsSolInt) {
    itsStepInSolInt = 0;
    ++itsBufferedSolInts;
    ++itsNSolInts;
  }

  if (itsBufferedSolInts == itsSolIntCount) {
    doSolve();

    // Clean up, prepare for next iteration
    itsAvgTime = 0;
    itsBufferedSolInts = 0;
    itsVisInInterval.assign(itsVisInInterval.size(),
                            std::pair<size_t, size_t>(0, 0));
    itsWeightsPerAntenna.assign(itsWeightsPerAntenna.size(), 0.0);

    for (size_t dir = 0; dir < itsResultSteps.size(); ++dir) {
      itsResultSteps[dir]->clear();
    }
    sol_ints_.clear();
  }

  ++itsTimeStep;
  itsTimer.stop();

  return false;
}

void DDECal::doPrepare(const DPBuffer& bufin, size_t sol_int, size_t step) {
  // UVW flagging happens on the copy of the buffer
  // These flags are later restored and therefore not written
  itsUVWFlagStep.process(bufin);

  itsTimerPredict.start();

  itsThreadPool->For(0, itsSteps.size(), [&](size_t dir, size_t) {
    itsSteps[dir]->process(bufin);
  });

  // Handle weights and flags
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nCr = 4;

  size_t nchanblocks = itsChanBlockFreqs.size();

  double weightFactor =
      1. / (nCh * (info().antennaUsed().size() - 1) * nCr * itsSolInt);

  for (size_t bl = 0; bl < nBl; ++bl) {
    size_t chanblock = 0, ant1 = info().antennaMap()[info().getAnt1()[bl]],
           ant2 = info().antennaMap()[info().getAnt2()[bl]];
    for (size_t ch = 0; ch < nCh; ++ch) {
      if (ch == itsChanBlockStart[chanblock + 1]) {
        chanblock++;
      }
      for (size_t cr = 0; cr < nCr; ++cr) {
        const size_t index = (bl * nCh + ch) * nCr + cr;
        itsVisInInterval[chanblock].second++;  // total nr of vis
        if (bufin.getFlags().data()[index]) {
          // Flagged points: set weight to 0
          sol_ints_[sol_int].WeightPtrs()[step][index] = 0;
        } else {
          // Add this weight to both involved antennas
          double weight = bufin.getWeights().data()[index];
          itsWeightsPerAntenna[ant1 * nchanblocks + chanblock] += weight;
          itsWeightsPerAntenna[ant2 * nchanblocks + chanblock] += weight;
          if (weight != 0.0) {
            itsVisInInterval[chanblock].first++;  // unflagged nr of vis
          }
        }
      }
    }
  }

  for (auto& weight : itsWeightsPerAntenna) {
    weight *= weightFactor;
  }

  itsTimerPredict.stop();

  itsAvgTime += itsAvgTime + bufin.getTime();
}

void DDECal::writeSolutions() {
  itsTimer.start();
  itsTimerWrite.start();

  unsigned int nSolTimes = (info().ntime() + itsSolInt - 1) / itsSolInt;
  unsigned int nDir = itsDirections.size();
  assert(nSolTimes == itsSols.size());
  std::vector<double> solTimes(nSolTimes);
  double starttime = info().startTime();
  for (unsigned int t = 0; t < nSolTimes; ++t) {
    solTimes[t] = starttime + (t + 0.5) * info().timeInterval() * itsSolInt;
  }

  if (itsConstraintSols[0].empty()) {
    // Record the actual iterands of the solver, not constraint results

    unsigned int nPol;

    std::vector<string> polarizations;
    if (itsMode == GainCal::DIAGONAL || itsMode == GainCal::DIAGONALPHASE ||
        itsMode == GainCal::DIAGONALAMPLITUDE) {
      nPol = 2;
      polarizations.emplace_back("XX");
      polarizations.emplace_back("YY");
    } else if (itsMode == GainCal::FULLJONES) {
      polarizations.emplace_back("XX");
      polarizations.emplace_back("XY");
      polarizations.emplace_back("YX");
      polarizations.emplace_back("YY");
      nPol = 4;
    } else {
      nPol = 1;
    }

    unsigned int nSolChan = itsSols[0].size();
    assert(nSolChan == itsChanBlockFreqs.size());

    size_t nSt = info().antennaUsed().size();
    std::vector<DComplex> sols(nSolChan * nSt * nSolTimes * nDir * nPol);
    size_t i = 0;

    // For nPol=1, loop over pol runs just once
    // For nPol=2, it runs over values 0 and 2 (picking diagonal elements from 4
    // pols) For nPol=4, it runs over 0, 1, 2, 3
    unsigned int polIncr = (itsMode == GainCal::FULLJONES ? 1 : 3);
    unsigned int maxPol = (nPol > 1 ? 4 : 1);
    // Put solutions in a contiguous piece of memory
    for (unsigned int time = 0; time < nSolTimes; ++time) {
      for (unsigned int chan = 0; chan < nSolChan; ++chan) {
        for (unsigned int ant = 0; ant < nSt; ++ant) {
          for (unsigned int dir = 0; dir < nDir; ++dir) {
            for (unsigned int pol = 0; pol < maxPol; pol += polIncr) {
              assert(!itsSols[time].empty());
              assert(!itsSols[time][chan].empty());
              assert(time < itsSols.size());
              assert(chan < itsSols[time].size());
              assert(ant * nDir * maxPol + dir * maxPol + pol <
                     itsSols[time][chan].size());
              assert(i < sols.size());
              sols[i] =
                  itsSols[time][chan][ant * nDir * maxPol + dir * maxPol + pol];
              ++i;
            }
          }
        }
      }
    }
    std::vector<schaapcommon::h5parm::AxisInfo> axes;
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("time", itsSols.size()));
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("freq", nSolChan));
    axes.emplace_back(
        schaapcommon::h5parm::AxisInfo("ant", info().antennaUsed().size()));
    axes.emplace_back(schaapcommon::h5parm::AxisInfo("dir", nDir));
    if (nPol > 1) {
      axes.emplace_back(schaapcommon::h5parm::AxisInfo("pol", nPol));
    }

    string historyString = "CREATE by DPPP\n" + DPPPVersion::AsString() + "\n" +
                           "step " + itsName + " in parset: \n" +
                           itsParsetString;
    unsigned int numsols = 1;
    // For [scalar]complexgain, store two soltabs: phase and amplitude
    if (itsMode == GainCal::DIAGONAL || itsMode == GainCal::SCALARCOMPLEXGAIN ||
        itsMode == GainCal::FULLJONES) {
      numsols = 2;
    }
    for (unsigned int solnum = 0; solnum < numsols; ++solnum) {
      string solTabName;
      schaapcommon::h5parm::SolTab soltab;
      switch (itsMode) {
        case GainCal::SCALARPHASE:
        case GainCal::DIAGONALPHASE:
        case GainCal::FULLJONES:
          if (solnum == 0) {
            solTabName = "phase000";
            soltab = itsH5Parm.CreateSolTab(solTabName, "phase", axes);
            soltab.SetComplexValues(sols, std::vector<double>(), false,
                                    historyString);
          } else {
            solTabName = "amplitude000";
            soltab = itsH5Parm.CreateSolTab(solTabName, "amplitude", axes);
            soltab.SetComplexValues(sols, std::vector<double>(), true,
                                    historyString);
          }
          break;
        case GainCal::SCALARCOMPLEXGAIN:
        case GainCal::DIAGONAL:
          if (solnum == 0) {
            solTabName = "phase000";
            soltab = itsH5Parm.CreateSolTab(solTabName, "phase", axes);
            soltab.SetComplexValues(sols, std::vector<double>(), false,
                                    historyString);
          } else {
            solTabName = "amplitude000";
            soltab = itsH5Parm.CreateSolTab(solTabName, "amplitude", axes);
            soltab.SetComplexValues(sols, std::vector<double>(), true,
                                    historyString);
          }
          break;
        case GainCal::SCALARAMPLITUDE:
        case GainCal::DIAGONALAMPLITUDE:
          solTabName = "amplitude000";
          soltab = itsH5Parm.CreateSolTab(solTabName, "amplitude", axes);
          soltab.SetComplexValues(sols, std::vector<double>(), true,
                                  historyString);
          break;
        default:
          throw std::runtime_error("Constraint should have produced output");
      }

      // Tell H5Parm which antennas were used
      std::vector<std::string> antennaUsedNames(info().antennaUsed().size());
      for (unsigned int i = 0; i < info().antennaUsed().size(); ++i) {
        antennaUsedNames[i] = info().antennaNames()[info().antennaUsed()[i]];
      }
      soltab.SetAntennas(antennaUsedNames);
      soltab.SetSources(getDirectionNames());

      if (nPol > 1) {
        soltab.SetPolarizations(polarizations);
      }

      soltab.SetFreqs(itsChanBlockFreqs);
      soltab.SetTimes(solTimes);
    }  // solnums loop
  } else {
    // Record the Constraint::Result in the H5Parm

    unsigned int nConstraints = itsConstraintSols[0].size();

    for (unsigned int constraintNum = 0; constraintNum < nConstraints;
         ++constraintNum) {
      // Number of solution names, e.g. 2 for "TEC" and "ScalarPhase"
      unsigned int nSolNames = itsConstraintSols[0][constraintNum].size();
      for (unsigned int solNameNum = 0; solNameNum < nSolNames; ++solNameNum) {
        // Get the result of the constraint solution at first time to get
        // metadata
        Constraint::Result firstResult =
            itsConstraintSols[0][constraintNum][solNameNum];

        std::vector<hsize_t> dims(firstResult.dims.size() +
                                  1);        // Add time dimension at beginning
        dims[0] = itsConstraintSols.size();  // Number of times
        size_t numSols = dims[0];
        for (unsigned int i = 1; i < dims.size(); ++i) {
          dims[i] = firstResult.dims[i - 1];
          numSols *= dims[i];
        }

        std::vector<string> firstaxesnames =
            StringUtil::tokenize(firstResult.axes, ",");

        std::vector<schaapcommon::h5parm::AxisInfo> axes;
        axes.emplace_back(
            schaapcommon::h5parm::AxisInfo("time", itsConstraintSols.size()));
        for (size_t axisNum = 0; axisNum < firstaxesnames.size(); ++axisNum) {
          axes.emplace_back(schaapcommon::h5parm::AxisInfo(
              firstaxesnames[axisNum], firstResult.dims[axisNum]));
        }

        // Put solutions in a contiguous piece of memory
        std::vector<double> sols(numSols);
        std::vector<double>::iterator nextpos = sols.begin();
        for (unsigned int time = 0; time < itsSols.size(); ++time) {
          if (itsConstraintSols[time].size() != itsConstraintSols[0].size())
            throw std::runtime_error(
                "Constraint " + std::to_string(constraintNum) +
                " did not produce a correct output at time step " +
                std::to_string(time) + ": got " +
                std::to_string(itsConstraintSols[time].size()) +
                " results, expecting " +
                std::to_string(itsConstraintSols[0].size()));
          nextpos = std::copy(
              itsConstraintSols[time][constraintNum][solNameNum].vals.begin(),
              itsConstraintSols[time][constraintNum][solNameNum].vals.end(),
              nextpos);
        }

        // Put solution weights in a contiguous piece of memory
        std::vector<double> weights;
        if (!itsConstraintSols[0][constraintNum][solNameNum].weights.empty()) {
          weights.resize(numSols);
          std::vector<double>::iterator nextpos = weights.begin();
          for (unsigned int time = 0; time < itsSols.size(); ++time) {
            nextpos =
                std::copy(itsConstraintSols[time][constraintNum][solNameNum]
                              .weights.begin(),
                          itsConstraintSols[time][constraintNum][solNameNum]
                              .weights.end(),
                          nextpos);
          }
        }

        string solTabName = firstResult.name + "000";
        schaapcommon::h5parm::SolTab soltab =
            itsH5Parm.CreateSolTab(solTabName, firstResult.name, axes);
        soltab.SetValues(sols, weights,
                         "CREATE by DPPP\n" + DPPPVersion::AsString() + "\n" +
                             "step " + itsName + " in parset: \n" +
                             itsParsetString);

        // Tell H5Parm which antennas were used
        std::vector<std::string> antennaUsedNames(info().antennaUsed().size());
        for (unsigned int i = 0; i < info().antennaUsed().size(); ++i) {
          antennaUsedNames[i] = info().antennaNames()[info().antennaUsed()[i]];
        }
        soltab.SetAntennas(antennaUsedNames);

        soltab.SetSources(getDirectionNames());

        if (soltab.HasAxis("pol")) {
          std::vector<string> polarizations;
          switch (soltab.GetAxis("pol").size) {
            case 2:
              polarizations.emplace_back("XX");
              polarizations.emplace_back("YY");
              break;
            case 4:
              polarizations.emplace_back("XX");
              polarizations.emplace_back("XY");
              polarizations.emplace_back("YX");
              polarizations.emplace_back("YY");
              break;
            default:
              throw std::runtime_error(
                  "No metadata for numpolarizations = " +
                  std::to_string(soltab.GetAxis("pol").size));
          }

          soltab.SetPolarizations(polarizations);
        }

        // Set channel to frequencies
        // Do not use itsChanBlockFreqs, because constraint may have changed
        // size
        unsigned int nChannelBlocks = 1;
        if (soltab.HasAxis("freq")) {
          nChannelBlocks = soltab.GetAxis("freq").size;
        }
        std::vector<double> chanBlockFreqs;

        chanBlockFreqs.resize(nChannelBlocks);
        for (size_t chBlock = 0; chBlock != nChannelBlocks; ++chBlock) {
          const size_t channelIndexStart =
                           chBlock * info().nchan() / nChannelBlocks,
                       channelIndexEnd =
                           (chBlock + 1) * info().nchan() / nChannelBlocks,
                       curChannelBlockSize =
                           channelIndexEnd - channelIndexStart;
          double meanfreq =
              std::accumulate(info().chanFreqs().data() + channelIndexStart,
                              info().chanFreqs().data() + channelIndexEnd,
                              0.0) /
              curChannelBlockSize;
          chanBlockFreqs[chBlock] = meanfreq;
        }

        soltab.SetFreqs(chanBlockFreqs);

        soltab.SetTimes(solTimes);
      }
    }
  }

  itsTimerWrite.stop();
  itsTimer.stop();
}

void DDECal::finish() {
  itsTimer.start();

  if (itsStepInSolInt > 0) {
    sol_ints_[itsBufferedSolInts].Fit();
  }

  if (sol_ints_.size() > 0) {
    doSolve();
  }

  if (!itsOnlyPredict) writeSolutions();

  sol_ints_.clear();
  itsTimer.stop();

  // Let the next steps finish.
  getNextStep()->finish();
}

void DDECal::storeModelData(
    const std::vector<std::vector<std::complex<float>*>>& input_model_data) {
  const size_t n_times = input_model_data.size();
  const size_t n_directions = itsSteps.size();
  const size_t n_visibilities =
      info().nbaselines() * info().nchan() * info().ncorr();
  // The model data might already have data in it. By using resize(),
  // reallocation between timesteps is avoided.
  itsModelData.resize(n_times);
  for (size_t timestep = 0; timestep != n_times; ++timestep) {
    itsModelData[timestep].resize(n_directions);
    for (size_t dir = 0; dir < n_directions; ++dir) {
      itsModelData[timestep][dir].assign(
          input_model_data[timestep][dir],
          input_model_data[timestep][dir] + n_visibilities);
    }
  }
}

void DDECal::subtractCorrectedModel(bool fullJones, size_t bufferIndex) {
  // The original data is still in the data buffers (the solver
  // doesn't change those). Here we apply the solutions to all the model data
  // directions and subtract them from the data.
  std::vector<std::vector<DComplex>>& solutions =
      itsSols[sol_ints_[bufferIndex].NSolution()];
  const size_t nBl = info().nbaselines();
  const size_t nCh = info().nchan();
  const size_t nDir = itsDirections.size();
  for (size_t time = 0; time != sol_ints_[bufferIndex].DataPtrs().size();
       ++time) {
    std::complex<float>* data = sol_ints_[bufferIndex].DataPtrs()[time];
    const std::vector<std::vector<std::complex<float>>>& modelData =
        itsModelData[time];
    for (size_t bl = 0; bl < nBl; ++bl) {
      const size_t ant1 = info().getAnt1()[bl];
      const size_t ant2 = info().getAnt2()[bl];
      size_t chanblock = 0;

      for (size_t ch = 0; ch < nCh; ++ch) {
        if (ch == itsChanBlockStart[chanblock + 1]) {
          chanblock++;
        }
        const size_t index = (bl * nCh + ch) * 4;
        if (itsOnlyPredict) {
          aocommon::MC2x2 value(aocommon::MC2x2::Zero());

          for (size_t dir = 0; dir != nDir; ++dir)
            value += aocommon::MC2x2(&modelData[dir][index]);

          for (size_t cr = 0; cr < 4; ++cr) data[index + cr] = value[cr];
        } else {
          aocommon::MC2x2 value(aocommon::MC2x2::Zero());
          for (size_t dir = 0; dir != nDir; ++dir) {
            if (fullJones) {
              aocommon::MC2x2 sol1(
                  &solutions[chanblock][(ant1 * nDir + dir) * 4]),
                  sol2(&solutions[chanblock][(ant2 * nDir + dir) * 4]);
              value += sol1.Multiply(aocommon::MC2x2(&modelData[dir][index]))
                           .MultiplyHerm(sol2);
            } else {
              std::complex<double> solfactor(
                  solutions[chanblock][ant1 * nDir + dir] *
                  std::conj(solutions[chanblock][ant2 * nDir + dir]));
              for (size_t cr = 0; cr < 4; ++cr)
                value[cr] += solfactor *
                             std::complex<double>(modelData[dir][index + cr]);
            }
          }
          for (size_t cr = 0; cr < 4; ++cr) data[index + cr] -= value[cr];
        }
      }  // channel loop
    }    // bl loop
  }      // time loop
}

}  // namespace DPPP
}  // namespace DP3
