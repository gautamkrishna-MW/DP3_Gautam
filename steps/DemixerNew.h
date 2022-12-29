// DemixerNew.h: DP3 step class to subtract A-team sources in adaptive way
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step class to subtract A-team sources in adaptive way
/// @author Ger van Diepen

#ifndef DP3_STEPS_DEMIXERNEW_H_
#define DP3_STEPS_DEMIXERNEW_H_

#include <ostream>

#include "Filter.h"

#include "../base/DemixInfo.h"
#include "../base/DemixWorker.h"
#include "../common/ParameterSet.h"
#include "../parmdb/ParmDB.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step class to subtract A-team sources in adaptive way

/// This class is a Step class to subtract the strong A-team sources
/// in a smart way (algorithm developed by Reinout van Weeren).
/// It operates as follows:
/// <ol>
/// <li> Per time window (default 2 min) demixing is done separately.
/// <li> Using a simple Ateam model (the A and other strong sources)
///      and the beam model, the StokesI flux is predicted per source.
///      The antennae of baselines with flux>threshold are counted.
///      Only an antenna counted in at least min_antenna baselines is
///      solved for. If not such antennae exist, the source is ignored.
/// <li> The target is predicted as well. Note that an Ateam source within
///      the target field is removed from the Ateam model and added/replaced
///      in the target model.
///      The target model is usually created using gsm.py.
/// <li> For the core baselines the ratio target/Ateam flux is used to
///      determine if the target has to be included, ignored or deprojected.
/// </ol>
/// It is based on the demixing.py script made by Bas vd Tol and operates
/// per time chunk as follows:
/// <ul>
///  <li> The data are phase-shifted and averaged for each source.
///  <li> Demixing is done using the combined results.
///  <li> For each source a BBS solve, smooth, and predict is done.
///  <li> The predicted results are subtracted from the averaged data.
/// </ul>

class DemixerNew : public Step {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  DemixerNew(const common::ParameterSet&, const string& prefix);

  common::Fields getRequiredFields() const override;

  common::Fields getProvidedFields() const override;

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the counts.
  void showCounts(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

 private:
  /// Process the data collected in itsBuf.
  void processData();

  /// Export the solutions to a ParmDB.
  void writeSolutions(double startTime, int ntime);

  /// Add the mean and M2 (square of differences) of a part in a
  /// numerically stable way.
  void addMeanM2(double& mean, double& m2, size_t& nr, double partmean,
                 double partm2, size_t partnr) const;

  /// Show a statistic.
  void showStat(std::ostream& os, double n, double ntot,
                const std::string& str1, const std::string& str2) const;

  /// Show a percentage with 1 decimal.
  void showPerc1(std::ostream& os, float perc) const;

  std::string itsName;
  base::DemixInfo itsDemixInfo;
  std::string itsInstrumentName;
  std::shared_ptr<parmdb::ParmDB> itsParmDB;
  Filter itsFilter;  ///< only used for getInfo()
  std::vector<base::DemixWorker> itsWorkers;
  std::vector<base::DPBuffer> itsBufIn;
  std::vector<base::DPBuffer> itsBufOut;
  std::vector<std::vector<double>>
      itsSolutions;                         ///< all solutions in a time window
  std::map<std::string, int> itsParmIdMap;  ///< -1 = new parm name
  unsigned int itsNTime;
  unsigned int itsNTimeOut;
  unsigned int itsNChunk;

  common::NSTimer itsTimer;
  common::NSTimer itsTimerDemix;
  common::NSTimer itsTimerDump;  ///< writeSolutions
  common::NSTimer itsTimerNext;  ///< next step (parallel to writeSolutions)
};

}  // namespace steps
}  // namespace dp3

#endif
