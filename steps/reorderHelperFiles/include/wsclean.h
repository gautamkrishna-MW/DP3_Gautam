#ifndef WSCLEAN_H
#define WSCLEAN_H

#include <aocommon/image.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>
#include <aocommon/imagecoordinates.h>

#include <schaapcommon/facets/facet.h>
#include <dp3/base/DPBuffer.h>

#include "imagingtable.h"
#include "observationinfo.h"
#include "outputchannelinfo.h"
#include "weightmode.h"
#include "partition.h"
#include "settings.h"

#include <optional>
#include <set>
#include <iostream>

class PrimaryBeam;

namespace schaapcommon {
namespace facets {
class FacetImage;
}
}  // namespace schaapcommon

class WSClean {
 public:
  WSClean();
  ~WSClean();

  Settings& GetSettings() { return _settings; }
  const Settings& GetSettings() const { return _settings; }
  void ResetSettings() { _settings = Settings(); }
  void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }

  void RunClean();

  /**
   * Entry point for performing a single prediction for an existing model image.
   *
   * In case of a facet-based prediction, the provided model images are assumed
   * to have the same size, so that the image size of the full image can be
   * inferred from the first entry in the _imagingTable in an early stage.
   */
  void RunPredict();
  void performReordering(bool isPredictMode);

  PartitionClass partitionObj;

 private:
  
  /**
   * Returns true when gridding is done with a-terms. This can either
   * be enabled by setting the gridWithBeam setting to true or by providing
   * an aterm config file. */
  bool griddingUsesATerms() const {
    return _settings.gridWithBeam || !_settings.atermConfigFilename.empty();
  }

  /**
   * True when the imaging uses any of the methods to apply a beam.
   * A beam can be applied through facetting (with solutions or beam),
   * through gridding with the beam using IDG or by correcting for the beam
   * in image space after imaging.
   */
  bool usesBeam() const {
    return _settings.applyPrimaryBeam || _settings.applyFacetBeam ||
           !_settings.facetSolutionFiles.empty() || griddingUsesATerms();
  }

  ObservationInfo getObservationInfo() const;
  /**
   * Add the phase shift of a facet
   * @param entry entry. If its facet is null, nothing happens.
   * @param l_shift is updated.
   * @param m_shift is updated.
   */

  std::pair<double, double> getLMShift() const;

  MSSelection selectInterval(MSSelection& fullSelection, size_t intervalIndex);

  void makeImagingTable(size_t outputIntervalIndex);
  void makeImagingTableEntry(const std::vector<aocommon::ChannelInfo>& channels,
                             size_t outIntervalIndex, size_t outChannelIndex,
                             ImagingTableEntry& entry);
  void makeImagingTableEntryChannelSettings(
      const std::vector<aocommon::ChannelInfo>& channels,
      size_t outIntervalIndex, size_t outChannelIndex, size_t nOutChannels,
      ImagingTableEntry& entry);
  void addPolarizationsToImagingTable(ImagingTableEntry& templateEntry);
  void addFacetsToImagingTable(ImagingTableEntry& templateEntry,
                               const size_t facet_count);
  void updateFacetsInImagingTable(
      const std::vector<std::shared_ptr<schaapcommon::facets::Facet>>& facets,
      bool updateDdPsfs);
  
  /**
   * Determines if IDG uses diagonal instrumental or full instrumental
   * polarizations.
   */
  aocommon::PolarizationEnum getProviderPolarization(
      aocommon::PolarizationEnum entry_polarization) const {
    if (_settings.gridderType == GridderType::IDG) {
      if (_settings.polarizations.size() == 1 &&
          *_settings.polarizations.begin() == aocommon::Polarization::StokesI) {
        if ((_settings.ddPsfGridWidth > 1 || _settings.ddPsfGridHeight > 1) &&
            _settings.gridWithBeam) {
          return aocommon::Polarization::StokesI;
        } else {
          return aocommon::Polarization::DiagonalInstrumental;
        }
      } else {
        return aocommon::Polarization::Instrumental;
      }
    } else if (_settings.diagonalSolutions) {
      return aocommon::Polarization::DiagonalInstrumental;
    } else {
      return entry_polarization;
    }
  }

  bool DataDescIdIsUsed(size_t ms_index, size_t data_desc_id) const {
    const size_t band_index = _msBands[ms_index].GetBandIndex(data_desc_id);
    // An empty selection means that all bands are selected
    return _settings.spectralWindows.empty() ||
           _settings.spectralWindows.find(band_index) !=
               _settings.spectralWindows.end();
  }

  MSSelection _globalSelection;
  std::string _commandLine;

  Settings _settings;

  std::vector<OutputChannelInfo> _infoPerChannel;
  OutputChannelInfo _infoForMFS;
  
  // std::vector<PartitionedMS::Handle> _partitionedMSHandles;
  std::vector<aocommon::MultiBandData> _msBands;
  ImagingTable _imagingTable;
  ObservationInfo _observationInfo;
  std::size_t _facetCount;  // 0 means facets are not used.
  std::size_t _ddPsfCount;  // 0 means dd-psfs are not used.
  /// These contain the user-requested image shift values converted from ra,dec
  /// to l,m units
  /// @{
  double _l_shift;
  double _m_shift;
  /// @}

};

#endif
