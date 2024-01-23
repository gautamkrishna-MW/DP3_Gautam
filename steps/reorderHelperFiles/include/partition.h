#include "settings.h"
#include "msselection.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <map>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <aocommon/logger.h>
#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <boost/filesystem/path.hpp>


using aocommon::Logger;

// We will create some efficiently packed structs to fetch data with 1 read.
// This will reduce the count of file-reads that are made.
// We do not apply this on the classes/structs themselves, because this will
//  reduce the data access performance.
#pragma pack(push, 1)
struct MetaRecordBuffer {
  double u;
  double v;
  double w;
  double time;
  uint16_t antenna1;
  uint16_t antenna2;
  uint16_t fieldId;
};
struct PartHeaderBuffer {
  uint64_t channelCount;
  uint64_t channelStart;
  uint32_t dataDescId;
  bool hasModel;
};
struct MetaHeaderBuffer {
  double startTime;
  uint64_t selectedRowCount;
  uint32_t filenameLength;
};
#pragma pack(pop)

struct MetaHeader {
    double startTime = 0.0;
    uint64_t selectedRowCount = 0;
    uint32_t filenameLength = 0;
    void Read(std::istream& str) {
      MetaHeaderBuffer metaHeaderBuffer;
      str.read(reinterpret_cast<char*>(&metaHeaderBuffer),
               sizeof(MetaHeaderBuffer));
      startTime = metaHeaderBuffer.startTime;
      selectedRowCount = metaHeaderBuffer.selectedRowCount;
      filenameLength = metaHeaderBuffer.filenameLength;
    }
    void Write(std::ostream& str) const {
      MetaHeaderBuffer metaHeaderBuffer{startTime, selectedRowCount,
                                        filenameLength};
      str.write(reinterpret_cast<const char*>(&metaHeaderBuffer),
                sizeof(MetaHeaderBuffer));
    }
    static constexpr size_t BINARY_SIZE =
        sizeof(startTime) + sizeof(selectedRowCount) + sizeof(filenameLength);
    static_assert(BINARY_SIZE == 20);
  } _metaHeader;

  struct MetaRecord {
    double u = 0.0, v = 0.0, w = 0.0, time = 0.0;
    uint16_t antenna1 = 0, antenna2 = 0, fieldId = 0;
    static constexpr size_t BINARY_SIZE =
        sizeof(double) * 4 + sizeof(uint16_t) * 3;
    static_assert(BINARY_SIZE == 38);
    void Read(std::istream& str) {
      MetaRecordBuffer metaRecordBuffer;
      str.read(reinterpret_cast<char*>(&metaRecordBuffer),
               sizeof(MetaRecordBuffer));

      u = metaRecordBuffer.u;
      v = metaRecordBuffer.v;
      w = metaRecordBuffer.w;
      time = metaRecordBuffer.time;
      antenna1 = metaRecordBuffer.antenna1;
      antenna2 = metaRecordBuffer.antenna2;
      fieldId = metaRecordBuffer.fieldId;
    }
    void Write(std::ostream& str) const {
      MetaRecordBuffer metaRecordBuffer{u,        v,        w,      time,
                                        antenna1, antenna2, fieldId};
      str.write(reinterpret_cast<const char*>(&metaRecordBuffer),
                sizeof(MetaRecordBuffer));
    }
  };

  struct PartHeader {
    uint64_t channelCount = 0;
    uint64_t channelStart = 0;
    uint32_t dataDescId = 0;
    bool hasModel = false;
    static constexpr size_t BINARY_SIZE = sizeof(channelCount) +
                                          sizeof(channelStart) +
                                          sizeof(dataDescId) + sizeof(hasModel);
    static_assert(BINARY_SIZE == 21);
    void Read(std::istream& str) {
      PartHeaderBuffer partHeaderBuffer;
      str.read(reinterpret_cast<char*>(&partHeaderBuffer),
               sizeof(PartHeaderBuffer));
      channelCount = partHeaderBuffer.channelCount;
      channelStart = partHeaderBuffer.channelStart;
      dataDescId = partHeaderBuffer.dataDescId;
      hasModel = partHeaderBuffer.hasModel;
    }
    void Write(std::ostream& str) const {
      PartHeaderBuffer partHeaderBuffer{channelCount, channelStart, dataDescId,
                                        hasModel};
      str.write(reinterpret_cast<const char*>(&partHeaderBuffer),
                sizeof(PartHeaderBuffer));
    }
  } _partHeader;

struct ChannelRange {
    int dataDescId;
    size_t start, end;
    bool operator<(const ChannelRange& rhs) const {
        if (dataDescId < rhs.dataDescId) return true;
        if (dataDescId > rhs.dataDescId) return false;
        if (start < rhs.start) return true;
        if (start > rhs.start) return false;
        return end < rhs.end;
    }
};

struct PartitionFiles {
  std::unique_ptr<std::ofstream> data;
  std::unique_ptr<std::ofstream> weight;
  std::unique_ptr<std::ofstream> model;
};

size_t GetMaxChannels(
    const std::vector<ChannelRange>& channel_ranges) {
  size_t max_channels = 0;
  for (const ChannelRange& range : channel_ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

std::string getFilenamePrefix(const std::string& msPathStr,
                                             const std::string& tempDir) {
  boost::filesystem::path prefixPath;
  if (tempDir.empty())
    prefixPath = msPathStr;
  else {
    std::string msPathCopy(msPathStr);
    while (!msPathCopy.empty() && *msPathCopy.rbegin() == '/')
      msPathCopy.resize(msPathCopy.size() - 1);
    boost::filesystem::path msPath(msPathCopy);
    prefixPath = boost::filesystem::path(tempDir) / msPath.filename();
  }
  std::string prefix(prefixPath.string());
  while (!prefix.empty() && *prefix.rbegin() == '/')
    prefix.resize(prefix.size() - 1);
  return prefix;
}

std::string getPartPrefix(const std::string& msPathStr,
                                         size_t partIndex,
                                         aocommon::PolarizationEnum pol,
                                         size_t dataDescId,
                                         const std::string& tempDir) {
  std::string prefix = getFilenamePrefix(msPathStr, tempDir);

  std::ostringstream partPrefix;
  partPrefix << prefix << "-part";
  if (partIndex < 1000) partPrefix << '0';
  if (partIndex < 100) partPrefix << '0';
  if (partIndex < 10) partPrefix << '0';
  partPrefix << partIndex;
  partPrefix << "-";
  partPrefix << aocommon::Polarization::TypeToShortString(pol);
  partPrefix << "-b" << dataDescId;
  return partPrefix.str();
}

std::map<size_t, size_t> getDataDescIdMap(const std::vector<ChannelRange>& channels) {
  std::map<size_t, size_t> dataDescIds;
  size_t spwIndex = 0;
  for (const ChannelRange& range : channels) {
    if (dataDescIds.count(range.dataDescId) == 0) {
      dataDescIds.emplace(range.dataDescId, spwIndex);
      ++spwIndex;
    }
  }
  return dataDescIds;
}

string getMetaFilename(const string& msPathStr,
                                      const std::string& tempDir,
                                      size_t dataDescId) {
  std::string prefix = getFilenamePrefix(msPathStr, tempDir);

  std::ostringstream s;
  s << prefix << "-spw" << dataDescId << "-parted-meta.tmp";
  return s.str();
}


std::vector<aocommon::PolarizationEnum> GetMSPolarizations(size_t dataDescId, const casacore::MeasurementSet& ms) {
  // First get the polarization index corresponding with the data desc id
  casacore::MSDataDescription dataDescriptionTable = ms.dataDescription();
  casacore::ScalarColumn<int> polarizationIndexColumn(
      dataDescriptionTable, casacore::MSDataDescription::columnName(
                                casacore::MSDataDescription::POLARIZATION_ID));
  const size_t polarizationIndex = polarizationIndexColumn(dataDescId);
  casacore::MSPolarization polTable = ms.polarization();
  std::vector<aocommon::PolarizationEnum> pols;
  casacore::ArrayColumn<int> corrTypeColumn(
      polTable, casacore::MSPolarization::columnName(
                    casacore::MSPolarizationEnums::CORR_TYPE));

  // Now get the information corresponding with the polarization index
  casacore::Array<int> corrTypeVec(corrTypeColumn(polarizationIndex));
  for (casacore::Array<int>::const_contiter p = corrTypeVec.cbegin();
       p != corrTypeVec.cend(); ++p) {
    pols.push_back(aocommon::Polarization::AipsIndexToEnum(*p));
  }

  return pols;
}

std::map<size_t, std::vector<aocommon::PolarizationEnum>>
GetMSPolarizationsPerDataDescId(const std::vector<ChannelRange>& ranges,
    casacore::MeasurementSet& ms) {
  std::map<size_t, std::vector<aocommon::PolarizationEnum>>
      msPolarizationsPerDataDescId;
  for (const ChannelRange& range : ranges) {
    msPolarizationsPerDataDescId.emplace(
        range.dataDescId,
        GetMSPolarizations(range.dataDescId, ms));
  }
  return msPolarizationsPerDataDescId;
}


void PartitionMain(const string& msPath, const std::vector<ChannelRange>& channels,
    MSSelection& selection, const string& dataColumnName, bool includeModel,
    bool initialModelRequired, const Settings& settings) {
  const bool modelUpdateRequired = settings.modelUpdateRequired;
  std::set<aocommon::PolarizationEnum> polsOut;
  if (settings.gridderType == GridderType::IDG) {
    if (settings.polarizations.size() == 1) {
      if ((settings.ddPsfGridWidth > 1 || settings.ddPsfGridHeight > 1) &&
          settings.gridWithBeam) {
        polsOut.insert(aocommon::Polarization::StokesI);
      } else {
        polsOut.insert(aocommon::Polarization::DiagonalInstrumental);
      }
    } else {
      polsOut.insert(aocommon::Polarization::Instrumental);
    }
  } else if (settings.diagonalSolutions) {
    polsOut.insert(aocommon::Polarization::DiagonalInstrumental);
  } else {
    polsOut = settings.polarizations;
  }
  const size_t polarizationsPerFile =
      aocommon::Polarization::GetVisibilityCount(*polsOut.begin());
  const std::string& temporaryDirectory = settings.temporaryDirectory;

  const size_t channelParts = channels.size();

  if (channelParts != 1) {
    Logger::Debug << "Partitioning in " << channels.size() << " channels:";
    for (size_t i = 0; i != channels.size(); ++i)
      Logger::Debug << ' ' << channels[i].dataDescId << ':' << channels[i].start
                    << '-' << channels[i].end;
  }
  Logger::Debug << '\n';

  // Ordered as files[pol x channelpart]
  std::vector<PartitionFiles> files(channelParts * polsOut.size());

  const size_t max_channels = GetMaxChannels(channels);

  // Each data desc id needs a separate meta file because they can have
  // different uvws and other info.
  size_t fileIndex = 0;
  for (size_t part = 0; part != channelParts; ++part) {
    for (aocommon::PolarizationEnum p : polsOut) {
      PartitionFiles& f = files[fileIndex];
      std::string partPrefix = getPartPrefix(
          msPath, part, p, channels[part].dataDescId, temporaryDirectory);
      f.data = std::make_unique<std::ofstream>(partPrefix + ".tmp");
      f.weight = std::make_unique<std::ofstream>(partPrefix + "-w.tmp");
      if (initialModelRequired)
        f.model = std::make_unique<std::ofstream>(partPrefix + "-m.tmp");
      f.data->seekp(PartHeader::BINARY_SIZE, std::ios::beg);

      ++fileIndex;
    }
  }

  // This maps dataDescId to spw index.
  const std::map<size_t, size_t> selectedDataDescIds =
      getDataDescIdMap(channels);

//   std::unique_ptr<MsRowProviderBase> rowProvider;
//   if (settings.baselineDependentAveragingInWavelengths == 0.0) {
//     if (settings.simulateNoise) {
//       std::unique_ptr<NoiseMSRowProvider> noiseRowProvider(
//           new NoiseMSRowProvider(msPath, selection, selectedDataDescIds,
//                                  dataColumnName, initialModelRequired));
//       if (settings.simulatedBaselineNoiseFilename.empty())
//         noiseRowProvider->SetNoiseLevel(settings.simulatedNoiseStdDev);
//       else
//         noiseRowProvider->SetNoiseBaselineFile(
//             settings.simulatedBaselineNoiseFilename);
//       rowProvider = std::move(noiseRowProvider);
//     } else
//       rowProvider = MakeMsRowProvider(msPath, selection, selectedDataDescIds,
//                                       dataColumnName, initialModelRequired);
//   } else {
//     if (initialModelRequired)
//       throw std::runtime_error(
//           "Baseline-dependent averaging is enabled together with a mode that "
//           "requires the model data (e.g. -continue or -subtract-model). This "
//           "is not possible.");
//     rowProvider = std::make_unique<AveragingMSRowProvider>(
//         settings.baselineDependentAveragingInWavelengths, msPath, selection,
//         selectedDataDescIds, settings.fieldIds[0], dataColumnName,
//         initialModelRequired);
//   }
  casacore::MeasurementSet msObj(msPath);
  const std::map<size_t, std::vector<aocommon::PolarizationEnum>>
      msPolarizationsPerDataDescId =
          GetMSPolarizationsPerDataDescId(channels, msObj);
  const size_t nAntennas = msObj.antenna().nrow(); // Error Check
  const aocommon::MultiBandData bands(msObj);

  if (settings.parallelReordering == 1)
    Logger::Info << "Reordering " << msPath << " into " << channelParts << " x "
                 << polsOut.size() << " parts.\n";

  // Write header of meta file, one meta file for each data desc id
  // TODO rather than writing we can just skip and write later
  std::vector<std::unique_ptr<std::ofstream>> metaFiles(
      selectedDataDescIds.size());
  for (const std::pair<const size_t, size_t>& p : selectedDataDescIds) {
    const size_t dataDescId = p.first;
    const size_t spwIndex = p.second;
    std::string metaFilename =
        getMetaFilename(msPath, temporaryDirectory, dataDescId);
    metaFiles[spwIndex] = std::make_unique<std::ofstream>(metaFilename);
    MetaHeader metaHeader;
    metaHeader.selectedRowCount = 0;  // not yet known
    metaHeader.filenameLength = msPath.size();
    metaHeader.startTime = 0;//rowProvider->StartTime();
    metaHeader.Write(*metaFiles[spwIndex]);
    metaFiles[spwIndex]->write(msPath.c_str(), msPath.size());
    if (!metaFiles[spwIndex]->good())
      throw std::runtime_error("Error writing to temporary file " +
                               metaFilename);
  }

  // Write actual data
  std::vector<std::complex<float>> dataBuffer(polarizationsPerFile *
                                              max_channels);
  std::vector<float> weightBuffer(polarizationsPerFile * max_channels);

  casacore::Array<std::complex<float>> dataArray;
  casacore::Array<std::complex<float>> modelArray;
  casacore::Array<float> weightSpectrumArray;
  casacore::Array<bool> flagArray;

  size_t selectedRowsTotal = 0;
  aocommon::UVector<size_t> selectedRowCountPerSpwIndex(
      selectedDataDescIds.size(), 0);
//   while (!rowProvider->AtEnd()) {
  int i=0;
  while (i<100) {
    i++;
    
    MetaRecord meta;

    double time;
    uint32_t dataDescId, antenna1, antenna2, fieldId;

    // Error
    // rowProvider->ReadData(dataArray, flagArray, weightSpectrumArray, meta.u,
    //                       meta.v, meta.w, dataDescId, antenna1, antenna2,
    //                       fieldId, time);
    meta.antenna1 = antenna1;
    meta.antenna2 = antenna2;
    meta.fieldId = fieldId;
    meta.time = time;
    const size_t spwIndex = selectedDataDescIds.find(dataDescId)->second;
    ++selectedRowCountPerSpwIndex[spwIndex];
    ++selectedRowsTotal;
    std::ofstream& metaFile = *metaFiles[spwIndex];
    meta.Write(metaFile);
    if (!metaFile.good())
      throw std::runtime_error("Error writing to temporary file");

    // if (initialModelRequired) rowProvider->ReadModel(modelArray);

    fileIndex = 0;
    for (size_t part = 0; part != channelParts; ++part) {
      if (channels[part].dataDescId == int(dataDescId)) {
        const size_t partStartCh = channels[part].start;
        const size_t partEndCh = channels[part].end;
        const std::vector<aocommon::PolarizationEnum>& msPolarizations =
            msPolarizationsPerDataDescId.find(dataDescId)->second;

        for (aocommon::PolarizationEnum p : polsOut) {
          PartitionFiles& f = files[fileIndex];
        //   CopyData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations,
        //            dataArray, p);
          f.data->write(reinterpret_cast<char*>(dataBuffer.data()),
                        (partEndCh - partStartCh) *
                            sizeof(std::complex<float>) * polarizationsPerFile);
          if (!f.data->good())
            throw std::runtime_error("Error writing to temporary data file");

          if (initialModelRequired) {
            // CopyData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations,
            //          modelArray, p);
            f.model->write(reinterpret_cast<char*>(dataBuffer.data()),
                           (partEndCh - partStartCh) *
                               sizeof(std::complex<float>) *
                               polarizationsPerFile);
            if (!f.model->good())
              throw std::runtime_error(
                  "Error writing to temporary model data file");
          }

        //   CopyWeights(weightBuffer.data(), partStartCh, partEndCh,
        //               msPolarizations, dataArray, weightSpectrumArray,
        //               flagArray, p);
          f.weight->write(
              reinterpret_cast<char*>(weightBuffer.data()),
              (partEndCh - partStartCh) * sizeof(float) * polarizationsPerFile);
          if (!f.weight->good())
            throw std::runtime_error("Error writing to temporary weights file");
          ++fileIndex;
        }
      } else {
        fileIndex += polsOut.size();
      }
    }

    // rowProvider->NextRow();
  }
  Logger::Debug << "Total selected rows: " << selectedRowsTotal << '\n';
//   rowProvider->OutputStatistics();

  // Rewrite meta headers to include selected row count
  for (const std::pair<const size_t, size_t>& p : selectedDataDescIds) {
    const size_t spwIndex = p.second;
    MetaHeader metaHeader;
    metaHeader.selectedRowCount = selectedRowCountPerSpwIndex[spwIndex];
    metaHeader.filenameLength = msPath.size();
    metaHeader.startTime = 0;//rowProvider->StartTime();
    metaFiles[spwIndex]->seekp(0);
    metaHeader.Write(*metaFiles[spwIndex]);
    metaFiles[spwIndex]->write(msPath.c_str(), msPath.size());
  }

  // Write header to parts and write empty model files (if requested)
  PartHeader header;
  header.hasModel = includeModel;
  fileIndex = 0;
  dataBuffer.assign(max_channels * polarizationsPerFile, 0.0);
  
//   if (includeModel && !initialModelRequired && settings.parallelReordering == 1)
//     progress2 =std::make_unique<ProgressBar>("Initializing model visibilities");
  for (size_t part = 0; part != channelParts; ++part) {
    header.channelStart = channels[part].start,
    header.channelCount = channels[part].end - header.channelStart;
    header.dataDescId = channels[part].dataDescId;
    for (std::set<aocommon::PolarizationEnum>::const_iterator p =
             polsOut.begin();
         p != polsOut.end(); ++p) {
      PartitionFiles& f = files[fileIndex];
      f.data->seekp(0, std::ios::beg);
      header.Write(*f.data);
      if (!f.data->good())
        throw std::runtime_error("Error writing to temporary data file");

      f.data.reset();
      f.weight.reset();
      f.model.reset();
      ++fileIndex;

      // If model is requested, fill model file with zeros
      if (includeModel && !initialModelRequired) {
        std::string partPrefix = getPartPrefix(
            msPath, part, *p, header.dataDescId, temporaryDirectory);
        std::ofstream modelFile(partPrefix + "-m.tmp");
        const size_t selectedRowCount = selectedRowCountPerSpwIndex
            [selectedDataDescIds.find(channels[part].dataDescId)->second];
        for (size_t i = 0; i != selectedRowCount; ++i) {
          modelFile.write(reinterpret_cast<char*>(dataBuffer.data()),
                          header.channelCount * sizeof(std::complex<float>) *
                              polarizationsPerFile);
        //   if (progress2)
        //     progress2->SetProgress(part * selectedRowCount + i,
        //                            channelParts * selectedRowCount);
        }
      }
    }
  }
//   progress2.reset();
}