// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step reading from an MS
/// @author Ger van Diepen

#ifndef DP3_STEPS_MSREADER_H_
#define DP3_STEPS_MSREADER_H_

#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableIter.h>

#include <dp3/base/DPBuffer.h>

#include "InputStep.h"

#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"

namespace dp3 {
namespace steps {
/// @brief DP3 step reading from an MS

/// This class is a InputStep step reading the data from a MeasurementSet.
/// At the beginning it finds out the shape of the data; i.e., the
/// number of correlations, channels, baselines, and time slots.
/// It requires the data to be regularly shaped.
///
/// The object is constructed from the 'msin' keywords in the parset file.
/// Currently the following can be given:
/// <ul>
///  <li> msin: name of the MS
///  <li> msin.autoweight: calculate weights from autocorrelations? [no]
///  <li> msin.startchan: first channel to use [0]
///  <li> msin.nchan: number of channels to use [all]
///  <li> msin.useflag: use the existing flags? [yes]
///  <li> msin.datacolumn: the data column to use [DATA]
///  <li> msin.weightcolumn: the weights column to use [WEIGHT_SPECTRUM or
///           WEIGHT]
///  <li> msin.starttime: first time to use [first time in MS]
///  <li> msin.endtime: last time to use [last time in MS]
/// </ul>
///
/// If a time slot is missing, it is inserted with flagged data set to zero.
/// Missing time slots can also be detected at the beginning or end of the
/// MS by giving the correct starttime and endtime.
/// The correct UVW coordinates are calculated for inserted time slots.
///
/// The process function reads the fields needed by the steps' chain which can
/// be obtained calling getFieldsToRead(). In this way only the necessary data
/// is kept in memory.
///
/// The data columns are handled in the following way:
/// <table>
///  <tr>
///   <td>TIME</td>
///   <td>The time slot center of the current data (in MJD seconds).
///       It is assumed that all data have the same interval, which is
///       used to find missing time slots.
///   </td>
///  </tr>
///  <tr>
///   <td>DATA</td>
///   <td>The visibility data as [ncorr,nchan,nbaseline]. Only the
///       part given by startchan and nchan is read. If a time slot is
///       inserted, all its data are zero.
///   </td>
///  </tr>
///  <tr>
///   <td>FLAG</td>
///   <td>The data flags as [ncorr,nchan,nbaseline] (True is bad).
///       They are read from the FLAG column. If a FLAG_ROW is set, all
///       flags for that baseline will be set. Also the flag of data
///       containing NaN or infinite numbers will be set.
///       All flags of an inserted time slot are set.
///   </td>
///  </tr>
///  <tr>
///   <td>WEIGHT</td>
///   <td>The data weights as [ncorr,nchan,nbaseline]. Column
///       WEIGHT_SPECTRUM is used if present and containing valid data,
///       otherwise column WEIGHT is used. The weights of an inserted
///       time slot are set to 0.
///       If autoweight is on, the autocorrelations are used to
///       calculate proper weights.
///   </td>
///  </tr>
///  <tr>
///   <td>UVW</td>
///   <td>The UVW coordinates in meters as [3,nbaseline].
///       They are calculated for a missing time slot.
///   </td>
///  </tr>
///  <tr>
///   <td>FULLRESFLAG</td>
///   <td>For each baseline the LOFAR_FULL_RES_FLAG column is stored as
///       a uChar array with shape [orignchan/8, ntimeavg]. The bits
///       represent the flags. They are converted to a Bool array with shape
///       [orignchan, ntimeavg, nbaseline].
///       If column LOFAR_FULL_RES_FLAG is not present, the flags are used
///       and it is assumed that no averaging was done yet (thus ntimeavg=1
///       and orignchan=nchan).
///   </td>
///  </tr>
/// </table>

class MSReader : public InputStep {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  /// The missingData argument is for MultiMSReader.
  MSReader(const casacore::MeasurementSet& ms, const common::ParameterSet&,
           const std::string& prefix, bool missingData = false);

  /// Override InputStep::getProvidedFields, since autoweight settings may add a
  /// field.
  common::Fields getProvidedFields() const override;

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Add some data to the MeasurementSet written/updated.
  /// Do nothing.
  void addToMS(const string&) override{};

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// If needed, show the flag counts.
  void showCounts(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Read the UVW at the given row numbers into the buffer.
  void getUVW(const casacore::RefRows& rowNrs, double time, base::DPBuffer&);

  /// Read the weights at the given row numbers into the buffer.
  /// Note: the buffer must contain DATA if autoweighting is in effect.
  void getWeights(const casacore::RefRows& rowNrs, base::DPBuffer&);

  /// Read the fullRes flags (LOFAR_FULL_RES_FLAG) at the given row numbers
  /// into the buffer.
  /// If there is no such column, the flags are set to false and false is
  /// returned.
  bool getFullResFlags(const casacore::RefRows& rowNrs, base::DPBuffer&);

  /// Fill the fullResFlags field in the DPBuffer.
  /// If the LOFAR_FULL_RES_FLAG column is present in the MS
  /// - the full resolution flags are read from the input.
  /// Else:
  /// - the full resolution flags are filled using the flags in the
  /// buffer.
  void FillFullResFlags(base::DPBuffer& buf);

  /// Get the main MS table.
  const casacore::Table& table() const override { return itsMS; }

  /// Get the slicer in the FLAG and DATA column.
  const casacore::Slicer& colSlicer() const { return itsColSlicer; }

  /// Get the rownrs for meta info of missing time slots.
  /// It uses the rows of the first time slot.
  const casacore::Vector<common::rownr_t>& getBaseRowNrs() const {
    return itsBaseRowNrs;
  }

  /// Get the name of the MS.
  std::string msName() const override;

  /// Get the time information.
  double firstTime() const override { return itsFirstTime; }
  double lastTime() const override { return itsLastTime; }

  /// Get the baseline selection.
  const string& baselineSelection() const { return itsSelBL; }

  /// Is the data column missing?
  bool missingData() const { return itsMissingData; }

  /// Get the start channel.
  unsigned int startChan() const { return itsStartChan; }

  /// Tell if the input MS has LOFAR_FULL_RES_FLAG.
  bool hasFullResFlags() const { return itsHasFullResFlags; }

  /// Get access to the buffer.
  const base::DPBuffer& getBuffer() const { return itsBuffer; }

  /// Flags inf and NaN
  static void flagInfNaN(const casacore::Cube<casacore::Complex>& dataCube,
                         casacore::Cube<bool>& flagsCube,
                         base::FlagCounter& flagCounter);

 protected:
  /// Default constructor.
  MSReader() = default;

 private:
  /// Prepare the access to the MS.
  /// Return the first and last time and the interval.
  void prepare(double& firstTime, double& lastTime, double& interval);

  /// Do the rest of the preparation.
  void prepare2(int spectralWindow);

  /// Skip the first times in the MS in case a start time was given.
  /// If needed, it sets itsFirstTime properly.
  void skipFirstTimes();

  /// Calculate the UVWs for a missing time slot.
  void calcUVW(double time, base::DPBuffer&);

  /// Calculate the weights from the autocorrelations.
  void autoWeight(casacore::Cube<float>& weights, const base::DPBuffer& buf);

 protected:
  casacore::MeasurementSet itsMS;
  casacore::Table itsSelMS;  ///< possible selection of spw, baseline
  casacore::TableIterator itsIter;
  std::string itsDataColName{"DATA"};
  std::string itsFlagColName{"FLAG"};
  std::string itsWeightColName{"WEIGHT_SPECTRUM"};
  std::string itsModelColName{"MODEL_DATA"};
  std::string itsStartChanStr{"0"};  ///< startchan expression
  std::string itsNrChanStr{"0"};     ///< nchan expression
  std::string itsSelBL{};            ///< Baseline selection string
  bool itsNeedSort{false};           ///< sort needed on time,baseline?
  bool itsAutoWeight{false};         ///< calculate weights from autocorr?
  bool itsAutoWeightForce{false};    ///< always calculate weights?
  bool itsHasWeightSpectrum;
  bool itsUseFlags{true};
  bool itsUseAllChan{false};   ///< all channels (i.e. no slicer)?
  bool itsMissingData{false};  ///< allow missing data column?
  unsigned int itsNrBl{0};
  unsigned int itsNrCorr{0};
  unsigned int itsNrChan{0};
  unsigned int itsStartChan{0};
  /// tolerance for time comparison
  ///
  /// Can be negative to insert flagged time slots before start.
  double itsTimeTolerance{1e-2};
  double itsTimeInterval{0.0};
  double itsFirstTime{0.0};
  double itsLastTime{0.0};
  double itsNextTime{0.0};
  double itsLastMSTime{0.0};
  unsigned int itsNrRead{0};      ///< nr of time slots read from MS
  unsigned int itsNrInserted{0};  ///< nr of inserted time slots
  casacore::Slicer itsColSlicer;  ///< slice in corr,chan column
  casacore::Slicer itsArrSlicer;  ///< slice in corr,chan,bl array
  bool itsHasFullResFlags{false};
  base::DPBuffer itsBuffer;
  std::unique_ptr<base::UVWCalculator> itsUVWCalc;
  casacore::Vector<common::rownr_t>
      itsBaseRowNrs;  ///< rownrs for meta of missing times
  base::FlagCounter itsFlagCounter;
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
