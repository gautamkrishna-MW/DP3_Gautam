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
  void RunClean();
  void RunPredict();
  void performReordering(bool isPredictMode);

  PartitionClass partitionObj;
  Settings _settings;

 private:
  ObservationInfo getObservationInfo() const;
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

  bool DataDescIdIsUsed(size_t ms_index, size_t data_desc_id) const {
    const size_t band_index = _msBands[ms_index].GetBandIndex(data_desc_id);
    // An empty selection means that all bands are selected
    return _settings.spectralWindows.empty() ||
           _settings.spectralWindows.find(band_index) !=
               _settings.spectralWindows.end();
  }

  MSSelection _globalSelection;
  std::vector<OutputChannelInfo> _infoPerChannel;
  std::vector<aocommon::MultiBandData> _msBands;
  ImagingTable _imagingTable;
  ObservationInfo _observationInfo;
  std::size_t _facetCount;  // 0 means facets are not used.
  std::size_t _ddPsfCount;  // 0 means dd-psfs are not used.
  double _l_shift;
  double _m_shift;
};

#endif
