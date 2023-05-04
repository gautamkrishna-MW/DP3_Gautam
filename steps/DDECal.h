// DDE.h: DP3 step class to calibrate direction dependent gains
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to apply a calibration correction to the data.
/// @author Tammo Jan Dijkema

#ifndef DP3_STEPS_DDECAL_H_
#define DP3_STEPS_DDECAL_H_

#include <fstream>
#include <string>
#include <vector>

#include <aocommon/threadpool.h>

#include "../common/ParameterSet.h"

#include "../ddecal/Settings.h"
#include "../ddecal/SolutionWriter.h"
#include "../ddecal/constraints/Constraint.h"
#include "../ddecal/gain_solvers/SolverBase.h"

#include "MultiResultStep.h"
#include "ResultStep.h"
#include "UVWFlagger.h"

namespace dp3 {
namespace steps {

/// @brief This class is a Step class to calibrate (direction independent)
/// gains.
class DDECal : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DDECal(const common::ParameterSet& parameterSet, const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  common::Fields getProvidedFields() const override {
    return (itsSettings.subtract || itsSettings.only_predict)
               ? kDataField
               : common::Fields();
  }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  void checkMinimumVisibilities(size_t bufferIndex);

  void flagChannelBlock(size_t cbIndex, size_t bufferIndex);

  /// Call the actual solver (called once per solution interval)
  void doSolve();

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  void updateInfo(const base::DPInfo&) override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

 private:
  void initializeColumnReaders(const common::ParameterSet&,
                               const string& prefix);
  void initializeIDG(const common::ParameterSet& parset, const string& prefix);
  void initializePredictSteps(const common::ParameterSet& parset,
                              const string& prefix);

  void setModelNextSteps(Step&, const std::string& direction,
                         const common::ParameterSet& parset,
                         const string& prefix) const;

  void doPrepare();

  /// Initializes solutions for a new solution interval.
  /// Based on progation settings, either copies the previous solution or
  /// writes default values to the new solution.
  /// @param buffer_index Index within the current solution interval set.
  void InitializeSolutions(size_t buffer_index);

  /// Write all solutions to an H5Parm file using itsSolutionWriter.
  void WriteSolutions();

  void subtractCorrectedModel(size_t bufferIndex);

  const ddecal::Settings itsSettings;

  /// The input data buffers for the current set of solution intervals.
  /// Maximum dimensions: itsSolIntCount x itsRequestedSolInt
  std::vector<std::vector<std::unique_ptr<base::DPBuffer>>> itsInputBuffers;
  /// Original flags of the input buffers for the current solution interval.
  /// This member is only used if itsUVWFlagger is active.
  /// Dimensions: ( solution_interval x step_within_interval x baseline x
  /// channel x correlation )
  xt::xtensor<bool, 5> itsOriginalFlags;

  /// The time of the current buffer (in case of solint, average time)
  double itsAvgTime;

  /// For each time, for each channel block, a vector of size nAntennas *
  /// SolverBase::NSolutions() * nPolarizations, with nPolarizations changing
  /// fastest.
  std::vector<std::vector<std::vector<casacore::DComplex>>> itsSols;
  std::vector<size_t> itsNIter,  // Number of iterations taken
      itsNApproxIter;

  /// For each time, for each constraint, a vector of results (e.g. tec and
  /// phase)
  std::vector<std::vector<std::vector<ddecal::Constraint::Result>>>
      itsConstraintSols;

  ddecal::SolutionWriter itsSolutionWriter;

  /// Number of timeslots to store per solution interval as requested
  /// by the user in the parset.
  size_t itsRequestedSolInt;
  /// For each direction, a number of solutions per solution interval
  std::vector<size_t> itsSolutionsPerDirection;
  size_t itsSolIntCount;  ///< Number of solution intervals to buffer
  /// Index of the first solution in the current solution interval set.
  size_t itsFirstSolutionIndex;
  size_t itsNChan;
  /// For each channel block, the nr of unflagged vis and the total nr of vis.
  std::vector<std::pair<size_t, size_t>> itsVisInInterval;
  /// For each channel block, the index in the channels at which this channel
  /// block starts.
  std::vector<size_t> itsChanBlockStart;
  std::vector<double> itsChanBlockFreqs;
  /// For each direction, a vector of patches.
  std::vector<std::vector<std::string>> itsDirections;
  /// Maps direction indices to the cluster central direction.
  std::vector<base::Direction> itsSourceDirections;

  /// First antenna for each baseline. Contains used antennas only.
  std::vector<int> itsAntennas1;
  /// Second antenna for each baseline. Contains used antennas only.
  std::vector<int> itsAntennas2;
  std::vector<double> itsWeightsPerAntenna;

  UVWFlagger itsUVWFlagStep;
  /// Result step for data after UV-flagging
  std::shared_ptr<ResultStep> itsDataResultStep;
  /// For each direction, the first step in the chain that computes the model.
  std::vector<std::shared_ptr<ModelDataStep>> itsSteps;
  /// For each direction, the required fields of the step chain.
  std::vector<common::Fields> itsRequiredFields;
  /// For each directions, a multiresultstep with all times.
  std::vector<std::shared_ptr<MultiResultStep>> itsResultSteps;

  /// Store the solution for later steps of processing in DPBuffer. Note: only
  /// works for 1 direction.
  bool itsStoreSolutionInBuffer;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerPredict;
  common::NSTimer itsTimerSolve;
  common::NSTimer itsTimerWrite;
  std::mutex itsMeasuresMutex;
  std::unique_ptr<ddecal::SolverBase> itsSolver;
  std::unique_ptr<aocommon::ThreadPool> itsThreadPool;
  std::unique_ptr<std::ofstream> itsStatStream;
};

}  // namespace steps
}  // namespace dp3

#endif
