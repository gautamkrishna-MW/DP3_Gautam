#include "../include/wsclean.h"
#include "../include/facetreader.h"
#include "../include/facetutil.h"
#include "../include/msselection.h"
#include "../include/progressbar.h"
#include "../include/imagingtable.h"

#include <aocommon/counting_semaphore.h>
#include <aocommon/dynamicfor.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/uvector.h>
#include <aocommon/units/angle.h>

#include <schaapcommon/facets/facetimage.h>
#include <schaapcommon/fft/resampler.h>
#include <schaapcommon/fft/restoreimage.h>
#include <schaapcommon/fitters/nlplfitter.h>

#include <algorithm>
#include <iostream>
#include <memory>

using aocommon::Image;
using aocommon::Logger;
using aocommon::Polarization;
using aocommon::PolarizationEnum;
using aocommon::units::Angle;

WSClean::WSClean()
    : _globalSelection(),
      _ddPsfCount(0) {}

WSClean::~WSClean() = default;

ObservationInfo WSClean::getObservationInfo() const {
  casacore::MeasurementSet ms(_settings.filenames[0]);
  ObservationInfo observationInfo =
      ReadObservationInfo(ms, _settings.fieldIds[0]);
  return observationInfo;
}

std::pair<double, double> WSClean::getLMShift() const {
  double l_shift = 0.0;
  double m_shift = 0.0;
  if (_settings.hasShift) {
    aocommon::ImageCoordinates::RaDecToLM(
        _settings.shiftRA, _settings.shiftDec, _observationInfo.phaseCentreRA,
        _observationInfo.phaseCentreDec, l_shift, m_shift);
  }
  return std::make_pair(l_shift, m_shift);
}

void WSClean::makeImagingTable(size_t outputIntervalIndex) {
  std::set<aocommon::ChannelInfo> channelSet;
  _msBands.assign(_settings.filenames.size(), aocommon::MultiBandData());
  for (size_t i = 0; i != _settings.filenames.size(); ++i) {
    casacore::MeasurementSet ms(_settings.filenames[i]);
    _msBands[i] = aocommon::MultiBandData(ms);
    std::set<size_t> dataDescIds = _msBands[i].GetUsedDataDescIds(ms);
    if (dataDescIds.size() != _msBands[i].DataDescCount()) {
      Logger::Debug << dataDescIds.size() << "/" << _msBands[i].DataDescCount()
                    << " spws are used of " << _settings.filenames[i] << '\n';
    }

    // Apply user selection: remove unselected spws
    if (!_settings.spectralWindows.empty()) {
      for (std::set<size_t>::iterator d = dataDescIds.begin();
           d != dataDescIds.end();) {
        if (_settings.spectralWindows.find(_msBands[i].GetBandIndex(*d)) ==
            _settings.spectralWindows.end())
          d = dataDescIds.erase(d);
        else
          ++d;
      }
    }
    // accumulate channel info
    for (const size_t dataDescId : dataDescIds) {
      bool increasing = true;
      if (_msBands[i][dataDescId].ChannelCount() >= 2) {
        increasing = _msBands[i][dataDescId].Channel(1) >
                     _msBands[i][dataDescId].Channel(0);
      }
      channelSet.insert(_msBands[i][dataDescId].Channel(0));
      for (size_t ch = 1; ch != _msBands[i][dataDescId].ChannelCount(); ++ch) {
        bool chanIncreasing = _msBands[i][dataDescId].Channel(ch) >
                              _msBands[i][dataDescId].Channel(ch - 1);
        if (chanIncreasing != increasing)
          throw std::runtime_error(
              "Your measurement set has an incorrect frequency axis: the "
              "channels do neither only increase nor only decrease in "
              "frequency");
        if (_msBands[i][dataDescId].Channel(ch) ==
            _msBands[i][dataDescId].Channel(ch - 1))
          throw std::runtime_error(
              "Your measurement set has an incorrect frequency axis: two "
              "adjacent channels had the same frequency. Channels should "
              "either strictly increase or strictly decrease in frequency.");
        channelSet.insert(_msBands[i][dataDescId].Channel(ch));
      }
    }
  }
  if (channelSet.size() < _settings.channelsOut) {
    std::ostringstream str;
    str << "Parameter '-channels-out' was set to an invalid value: "
        << _settings.channelsOut
        << " output channels requested, but combined in all specified "
           "measurement sets, there are only "
        << channelSet.size() << " unique channels.";
    throw std::runtime_error(str.str());
  }
  std::vector<aocommon::ChannelInfo> inputChannelFrequencies(channelSet.begin(),
                                                             channelSet.end());
  Logger::Debug << "Total nr of channels found in measurement sets: "
                << inputChannelFrequencies.size() << '\n';

  _imagingTable.Clear();

  ImagingTableEntry templateEntry;
  templateEntry.joinedGroupIndex = 0;
  templateEntry.squaredDeconvolutionIndex = 0;

  // for(size_t interval=0; interval!=_settings.intervalsOut; ++interval)
  //{
  for (size_t outChannelIndex = 0; outChannelIndex != _settings.channelsOut;
       ++outChannelIndex) {
    makeImagingTableEntry(inputChannelFrequencies, outputIntervalIndex,
                          outChannelIndex, templateEntry);
    templateEntry.outputChannelIndex = outChannelIndex;

    if (_settings.joinedFrequencyDeconvolution) {
      templateEntry.joinedGroupIndex = 0;
    }
    addPolarizationsToImagingTable(templateEntry);
  }
  //}
  _imagingTable.Update();
  _imagingTable.Print();
}

void WSClean::makeImagingTableEntry(
    const std::vector<aocommon::ChannelInfo>& channels, size_t outIntervalIndex,
    size_t outChannelIndex, ImagingTableEntry& entry) {
  size_t startCh, endCh;
  if (_settings.endChannel != 0) {
    if (_settings.endChannel > channels.size())
      throw std::runtime_error(
          "Bad channel selection -- more channels selected than available");
    startCh = _settings.startChannel;
    endCh = _settings.endChannel;
  } else {
    startCh = 0;
    endCh = channels.size();
  }
  std::vector<aocommon::ChannelInfo> groupChannels(channels.begin() + startCh,
                                                   channels.begin() + endCh);

  if (_settings.divideChannelFrequencies.empty()) {
    makeImagingTableEntryChannelSettings(groupChannels, outIntervalIndex,
                                         outChannelIndex, _settings.channelsOut,
                                         entry);
  } else {
    // We need to separately divide the channels into groups as specified and
    // call the freq division for the group corresponding with the
    // outChannelIndex.
    const size_t nSplits = _settings.divideChannelFrequencies.size();
    for (size_t i = 0; i != nSplits + 1; ++i) {
      const size_t outChannelStart = _settings.channelsOut * i / (nSplits + 1);
      const size_t outChannelEnd =
          _settings.channelsOut * (i + 1) / (nSplits + 1);
      if (outChannelIndex >= outChannelStart &&
          outChannelIndex < outChannelEnd) {
        double splitFreqLow =
            (i == 0) ? 0.0 : _settings.divideChannelFrequencies[i - 1];
        double splitFreqHigh = (i == nSplits)
                                   ? std::numeric_limits<double>::max()
                                   : _settings.divideChannelFrequencies[i];
        std::vector<aocommon::ChannelInfo> splittedChannels;
        for (const aocommon::ChannelInfo& channel : groupChannels) {
          if (channel.Frequency() >= splitFreqLow &&
              channel.Frequency() < splitFreqHigh)
            splittedChannels.emplace_back(channel);
        }
        size_t nOutChannels = outChannelEnd - outChannelStart;
        makeImagingTableEntryChannelSettings(splittedChannels, outIntervalIndex,
                                             outChannelIndex - outChannelStart,
                                             nOutChannels, entry);
      }
    }
  }

  if (_settings.spectralCorrection.empty())
    entry.siCorrection = 1.0;
  else {
    double bandwidthCentre =
        0.5 * (channels.front().Frequency() + channels.back().Frequency());
    double chCentralFrequency =
        0.5 * (entry.lowestFrequency + entry.highestFrequency);
    double chFlux = schaapcommon::fitters::NonLinearPowerLawFitter::Evaluate(
        chCentralFrequency, _settings.spectralCorrection,
        _settings.spectralCorrectionFrequency);
    double midFlux = schaapcommon::fitters::NonLinearPowerLawFitter::Evaluate(
        bandwidthCentre, _settings.spectralCorrection,
        _settings.spectralCorrectionFrequency);
    entry.siCorrection = midFlux / chFlux;
    if (outChannelIndex == 0)
      Logger::Debug << "SI correction for first channel: " << entry.siCorrection
                    << '\n';
    if (outChannelIndex + 1 == _settings.channelsOut)
      Logger::Debug << "SI correction for last channel: " << entry.siCorrection
                    << '\n';
  }

  entry.msData.resize(_settings.filenames.size());
  for (size_t msIndex = 0; msIndex != _settings.filenames.size(); ++msIndex) {
    entry.msData[msIndex].bands.resize(_msBands[msIndex].DataDescCount());
  }
}

void WSClean::performReordering(bool isPredictMode) {
  std::mutex mutex;

  partitionObj.includeModel = _settings.deconvolutionMGain != 1.0 || isPredictMode ||
                  _settings.subtractModel || _settings.continuedRun;
  partitionObj.initialModelRequired = _settings.subtractModel || _settings.continuedRun;

  if (_settings.parallelReordering != 1) Logger::Info << "Reordering...\n";

  aocommon::CountingSemaphore semaphore(_settings.parallelReordering);
  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, _settings.filenames.size(), [&](size_t msIndex) {
    aocommon::ScopedCountingSemaphoreLock semaphore_lock(semaphore);
    // std::vector<ChannelRange> channels;
    // The partIndex needs to increase per data desc ids and channel ranges
    std::map<PolarizationEnum, size_t> nextIndex;
    for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroupCount();
         ++sqIndex) {
      ImagingTable sqGroup = _imagingTable.GetSquaredGroup(sqIndex);
      ImagingTable::Groups facet_groups = sqGroup.FacetGroups(true);
      for (size_t fgIndex = 0; fgIndex != facet_groups.size(); ++fgIndex) {
        ImagingTable facetGroup = ImagingTable(facet_groups[fgIndex]);
        // The band information is determined from the first facet in the group.
        // After this, all facet entries inside the group are updated.
        const ImagingTableEntry& entry = facetGroup.Front();
        for (size_t d = 0; d != _msBands[msIndex].DataDescCount(); ++d) {
          MSSelection selection(_globalSelection);
          if (DataDescIdIsUsed(msIndex, d) &&
              selection.SelectMsChannels(_msBands[msIndex], d, entry)) {
            if (entry.polarization == *_settings.polarizations.begin()) {
              ChannelRange r;
              r.dataDescId = d;
              r.start = selection.ChannelRangeStart();
              r.end = selection.ChannelRangeEnd();
              partitionObj.channels.push_back(r);
            }
            for (ImagingTableEntry& facetEntry : facetGroup) {
              facetEntry.msData[msIndex].bands[d].partIndex =
                  nextIndex[entry.polarization];
            }
            ++nextIndex[entry.polarization];
          }
        }
      }
    }

    partitionObj.msPath = _settings.filenames[msIndex];
    partitionObj.preprocessPartition(_settings.dataColumnName, _settings);
    std::lock_guard<std::mutex> lock(mutex);
  });
}

void WSClean::makeImagingTableEntryChannelSettings(
    const std::vector<aocommon::ChannelInfo>& channels, size_t outIntervalIndex,
    size_t outChannelIndex, size_t nOutChannels, ImagingTableEntry& entry) {
  size_t chLowIndex, chHighIndex;
  if (_settings.divideChannelsByGaps) {
    std::multimap<double, size_t> gaps;
    for (size_t i = 1; i != channels.size(); ++i) {
      double left = channels[i - 1].Frequency();
      double right = channels[i].Frequency();
      gaps.emplace(right - left, i);
    }
    std::vector<size_t> orderedGaps;
    auto iter = gaps.rbegin();
    for (size_t i = 0; i != nOutChannels - 1; ++i) {
      if (iter == gaps.rend())
        throw std::runtime_error(
            "Channel gap division leads to invalid selection");
      orderedGaps.push_back(iter->second);
      ++iter;
    }
    std::sort(orderedGaps.begin(), orderedGaps.end());
    if (outChannelIndex == 0)
      chLowIndex = 0;
    else
      chLowIndex = orderedGaps[outChannelIndex - 1];
    if (outChannelIndex + 1 == nOutChannels)
      chHighIndex = channels.size() - 1;
    else
      chHighIndex = orderedGaps[outChannelIndex] - 1;
  } else {
    chLowIndex = outChannelIndex * channels.size() / nOutChannels;
    chHighIndex = (outChannelIndex + 1) * channels.size() / nOutChannels - 1;
    if (chLowIndex == chHighIndex + 1)
      throw std::runtime_error(
          "Too many output channels requested: output channel " +
          std::to_string(outChannelIndex) +
          " would be empty. Number of output channels requested: " +
          std::to_string(_settings.channelsOut) +
          ". Number of channels in the measurement set(s) available (after "
          "applying channel range selections and splits): " +
          std::to_string(channels.size()));
  }
  if (channels[chLowIndex].Frequency() > channels[chHighIndex].Frequency())
    std::swap(chLowIndex, chHighIndex);
  entry.inputChannelCount = chHighIndex + 1 - chLowIndex;
  entry.lowestFrequency = channels[chLowIndex].Frequency();
  entry.highestFrequency = channels[chHighIndex].Frequency();
  entry.bandStartFrequency =
      entry.lowestFrequency - channels[chLowIndex].Width() * 0.5;
  entry.bandEndFrequency =
      entry.highestFrequency + channels[chHighIndex].Width() * 0.5;
  entry.outputIntervalIndex = outIntervalIndex;
}

void WSClean::addPolarizationsToImagingTable(ImagingTableEntry& templateEntry) {
  for (PolarizationEnum p : _settings.polarizations) {
    const bool isFirstPol = (p == *_settings.polarizations.begin());
    templateEntry.polarization = p;
    if (p == Polarization::XY)
      templateEntry.imageCount = 2;
    else if (p == Polarization::YX)
      templateEntry.imageCount = 0;
    else
      templateEntry.imageCount = 1;

    if (_ddPsfCount && isFirstPol) {
      ImagingTableEntry ddPsfTemplateEntry(templateEntry);
      ddPsfTemplateEntry.isDdPsf = true;
      addFacetsToImagingTable(ddPsfTemplateEntry, _ddPsfCount);
    }
    addFacetsToImagingTable(templateEntry, _facetCount);

    if (!_settings.joinedPolarizationDeconvolution) {
      ++templateEntry.joinedGroupIndex;
      ++templateEntry.squaredDeconvolutionIndex;
    }
  }

  if (_settings.joinedPolarizationDeconvolution) {
    ++templateEntry.joinedGroupIndex;
    ++templateEntry.squaredDeconvolutionIndex;
  }
}

void WSClean::addFacetsToImagingTable(ImagingTableEntry& templateEntry,
                                      const size_t facet_count) {
  // Create a single entry (with facetIndex == 0) when facets are not used.
  const size_t facet_entry_count = std::max(facet_count, std::size_t(1));
  for (size_t f = 0; f != facet_entry_count; ++f) {
    auto entry = std::make_unique<ImagingTableEntry>(templateEntry);
    entry->facetIndex = f;
    entry->facet.reset();  // updateFacetsInImagingTable will set the facet.
    _imagingTable.AddEntry(std::move(entry));
  }
  ++templateEntry.facetGroupIndex;
}

void WSClean::updateFacetsInImagingTable(
    const std::vector<std::shared_ptr<schaapcommon::facets::Facet>>& facets,
    bool updateDdPsfs) {
  for (ImagingTableEntry& entry : _imagingTable) {
    if (entry.isDdPsf != updateDdPsfs) continue;
    assert(entry.facetIndex < facets.size());
    entry.facet = facets[entry.facetIndex];
    // Calculate phase center delta for entry
    entry.centreShiftX = entry.facet->GetUntrimmedBoundingBox().Centre().x -
                         _settings.trimmedImageWidth / 2;
    entry.centreShiftY = entry.facet->GetUntrimmedBoundingBox().Centre().y -
                         _settings.trimmedImageHeight / 2;
  }
}

MSSelection WSClean::selectInterval(MSSelection& fullSelection,
                                    size_t intervalIndex) {
  if (_settings.intervalsOut == 1)
    return fullSelection;
  else {
    size_t tS, tE;
    if (fullSelection.HasInterval()) {
      tS = fullSelection.IntervalStart();
      tE = fullSelection.IntervalEnd();
    } else {
      casacore::MeasurementSet ms(_settings.filenames[0]);
      Logger::Info << "Counting number of scans... ";
      Logger::Info.Flush();
      casacore::ScalarColumn<double> timeColumn(
          ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
      double time = timeColumn(0);
      size_t timestepIndex = 1;
      for (size_t row = 0; row != ms.nrow(); ++row) {
        if (time != timeColumn(row)) {
          ++timestepIndex;
          time = timeColumn(row);
        }
      }
      Logger::Info << "DONE (" << timestepIndex << ")\n";
      tS = 0;
      tE = timestepIndex;
      // Store the full interval in the selection, so that it doesn't need to
      // be determined again.
      fullSelection.SetInterval(tS, tE);
    }
    if (_settings.intervalsOut > tE - tS) {
      std::ostringstream str;
      str << "Invalid interval selection: " << _settings.intervalsOut
          << " intervals requested, but measurement set has only " << tE - tS
          << " intervals.";
      throw std::runtime_error(str.str());
    }
    MSSelection newSelection(fullSelection);
    newSelection.SetInterval(
        tS + (tE - tS) * intervalIndex / _settings.intervalsOut,
        tS + (tE - tS) * (intervalIndex + 1) / _settings.intervalsOut);
    return newSelection;
  }
}

void WSClean::RunPredict() {
  _observationInfo = getObservationInfo();
  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> facets;
  _facetCount = FacetReader::CountFacets(_settings.facetRegionFilename);
  std::tie(_l_shift, _m_shift) = getLMShift();
  makeImagingTable(0);
  _infoPerChannel.assign(_settings.channelsOut, OutputChannelInfo());
  _globalSelection = _settings.GetMSSelection();
  MSSelection fullSelection = _globalSelection;
  _globalSelection = selectInterval(fullSelection, 0);
}

void WSClean::RunClean() {
  _observationInfo = getObservationInfo();
  std::tie(_l_shift, _m_shift) = getLMShift();

  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> facets =
      FacetReader::ReadFacets(
          _settings.facetRegionFilename, _settings.trimmedImageWidth,
          _settings.trimmedImageHeight, _settings.pixelScaleX,
          _settings.pixelScaleY, _observationInfo.phaseCentreRA,
          _observationInfo.phaseCentreDec, _l_shift, _m_shift,
          _settings.imagePadding, _settings.gridderType == GridderType::IDG,
          _settings.GetFeatherSize());
  _facetCount = facets.size();

  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> dd_psfs;
  if (_settings.ddPsfGridWidth > 1 || _settings.ddPsfGridHeight > 1) {
    const schaapcommon::facets::Facet::InitializationData facet_data =
        CreateFacetInitializationData(
            _settings.trimmedImageWidth, _settings.trimmedImageHeight,
            _settings.pixelScaleX, _settings.pixelScaleY,
            _observationInfo.phaseCentreRA, _observationInfo.phaseCentreDec,
            _l_shift, _m_shift, _settings.imagePadding,
            _settings.gridderType == GridderType::IDG, 0);
    dd_psfs = CreateFacetGrid(facet_data, _settings.ddPsfGridWidth,
                              _settings.ddPsfGridHeight);
  }
  _ddPsfCount = dd_psfs.size();

  schaapcommon::facets::PixelPosition centerPixel(
      _settings.trimmedImageWidth / 2, _settings.trimmedImageHeight / 2);
  const bool hasCenter = std::any_of(
      facets.begin(), facets.end(),
      [&centerPixel](
          const std::shared_ptr<schaapcommon::facets::Facet>& facet) {
        // Point-in-poly test only evaluated if bounding box does
        // contain the centerPixel
        return facet->GetTrimmedBoundingBox().Contains(centerPixel) &&
               facet->Contains(centerPixel);
      });
  // FIXME: raise warning if facets do not cover the entire image, see AST-429

  // Center pixel should be present in one of the facets for the deconvolution
  if (!facets.empty() && _settings.deconvolutionIterationCount > 0 &&
      !hasCenter) {
    throw std::runtime_error(
        "The center pixel of the full image is not found in one of the facets. "
        "Make sure your facet file defines a facet that covers the center "
        "pixel of the main image.");
  }

  _globalSelection = _settings.GetMSSelection();
  MSSelection fullSelection = _globalSelection;

  makeImagingTable(0);
  if (!facets.empty()) updateFacetsInImagingTable(facets, false);
  if (!dd_psfs.empty()) updateFacetsInImagingTable(dd_psfs, true);

   _globalSelection = selectInterval(fullSelection, 0);
}
