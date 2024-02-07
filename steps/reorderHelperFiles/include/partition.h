
#ifndef PARTITION_H
#define PARTITION_H

#include <string>
#include <iostream>
#include "msselection.h"
#include "settings.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include <dp3/steps/Step.h>
#include <dp3/common/Types.h>

using dp3::base::DPInfo;
using dp3::steps::Step;

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

template <typename NumType>
  static bool IsCFinite(const std::complex<NumType>& c) {
    return std::isfinite(c.real()) && std::isfinite(c.imag());
  }

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

class PartitionClass
{
public:
  PartitionClass():
      msPath(),
      temporaryDirectory()
      {}
  ~PartitionClass(){}

  bool initialModelRequired;
  bool includeModel;
  std::string msPath;
  std::vector<ChannelRange> channels;
  dp3::base::DPInfo infoObj;

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
  };

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
      MetaRecordBuffer metaRecordBuffer{u, v, w, time, antenna1, antenna2, fieldId};
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
  };

  size_t GetMaxChannels(const std::vector<ChannelRange>& channel_ranges);
  string getMetaFilename(const string& msPathStr, const std::string& tempDir, size_t dataDescId);
  std::vector<aocommon::PolarizationEnum> GetMSPolarizations(size_t dataDescId, const casacore::MeasurementSet& ms);
  std::map<size_t, size_t> getDataDescIdMap(const std::vector<ChannelRange>& channels);

  std::string getFilenamePrefix(const std::string& msPathStr,
    const std::string& tempDir);

  std::string getPartPrefix(const std::string& msPathStr, size_t partIndex,
    aocommon::PolarizationEnum pol, size_t dataDescId, const std::string& tempDir);
  
  std::map<size_t, std::vector<aocommon::PolarizationEnum>>
  GetMSPolarizationsPerDataDescId(const std::vector<ChannelRange>& ranges, casacore::MeasurementSet& ms);

  void CopyData(std::complex<float>* dest, size_t startChannel,
    size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    aocommon::PolarizationEnum polOut);

  template <typename NumType>
  void CopyWeights(NumType* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut);

  void preprocessPartition(const string& dataColumnName, const Settings& settings);
  void processPartition(dp3::base::DPBuffer* buffer);
  void postprocess();

private:
  size_t polarizationsPerFile;
  size_t channelParts;
  size_t max_channels;
  std::string temporaryDirectory;
  std::vector<PartitionFiles> files;  
  std::vector<std::complex<float>> dataBuffer;
  std::set<aocommon::PolarizationEnum> polsOut;
  std::map<size_t, size_t> selectedDataDescIds;
  std::vector<std::unique_ptr<std::ofstream>> metaFiles;
  aocommon::UVector<size_t> selectedRowCountPerSpwIndex;
};


#endif