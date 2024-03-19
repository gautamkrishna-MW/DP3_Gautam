#include "../include/partition.h"

#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <random>

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

#include "../include/progressbar.h"

using aocommon::Logger;
using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::rownr_t;
using dp3::steps::Step;

PartitionClass::PartHeader _partHeader;
PartitionClass::MetaHeader _metaHeader;

PartitionClass::NoiseMap::NoiseMap(std::istream& stream) {
  size_t maxAnt = 0;
  std::string line;
  while (stream) {
    std::getline(stream, line);
    if (stream) {
      std::stringstream linestr(line);
      size_t ant1, ant2;
      linestr >> ant1 >> ant2;
      if (linestr) {
        float stddev;
        linestr >> stddev;
        if (!linestr) stddev = std::numeric_limits<float>::quiet_NaN();
        if (ant1 > ant2) std::swap(ant1, ant2);
        maxAnt = std::max(maxAnt, std::max(ant1, ant2));
        const bool isInserted =
            _map.emplace(std::make_pair(ant1, ant2), stddev).second;
        if (!isInserted)
          throw std::runtime_error(
              "Baseline " + std::to_string(ant1) + " x " +
              std::to_string(ant2) +
              " is specified twice in the noise baseline file");
      }
    }
  }
  aocommon::Logger::Info << "Read noise baseline file with " << _map.size()
                         << " rows and " << maxAnt + 1 << " antennas.\n";
}

float PartitionClass::NoiseMap::GetNoiseValue(size_t antenna1,
                                              size_t antenna2) const {
  size_t a1 = antenna1;
  size_t a2 = antenna2;
  if (a1 > a2) std::swap(a1, a2);
  auto iter = _map.find(std::make_pair(a1, a2));
  if (iter == _map.end())
    throw std::runtime_error(
        "The following baseline was not present in the baseline noise "
        "map: " +
        std::to_string(antenna1) + " x " + std::to_string(antenna2));
  return iter->second;
}

size_t PartitionClass::GetMaxChannels(
    const std::vector<ChannelRange>& channel_ranges) {
  size_t max_channels = 0;
  for (const ChannelRange& range : channel_ranges) {
    max_channels = std::max(max_channels, range.end - range.start);
  }
  return max_channels;
}

std::string PartitionClass::getFilenamePrefix(const std::string& msPathStr,
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

std::string PartitionClass::getPartPrefix(const std::string& msPathStr,
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

std::map<size_t, size_t> PartitionClass::getDataDescIdMap(
    const std::vector<ChannelRange>& channels) {
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

string PartitionClass::getMetaFilename(const string& msPathStr,
                                       const std::string& tempDir,
                                       size_t dataDescId) {
  std::string prefix = getFilenamePrefix(msPathStr, tempDir);

  std::ostringstream s;
  s << prefix << "-spw" << dataDescId << "-parted-meta.tmp";
  return s.str();
}

std::vector<aocommon::PolarizationEnum> PartitionClass::GetMSPolarizations(
    size_t dataDescId, const casacore::MeasurementSet& ms) {
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
PartitionClass::GetMSPolarizationsPerDataDescId(
    const std::vector<ChannelRange>& ranges, casacore::MeasurementSet& ms) {
  std::map<size_t, std::vector<aocommon::PolarizationEnum>>
      msPolarizationsPerDataDescId;
  for (const ChannelRange& range : ranges) {
    msPolarizationsPerDataDescId.emplace(
        range.dataDescId, GetMSPolarizations(range.dataDescId, ms));
  }
  return msPolarizationsPerDataDescId;
}

void PartitionClass::CopyData(
    std::complex<float>* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    aocommon::PolarizationEnum polOut) {
  const size_t polCount = polsIn.size();
  casacore::Array<std::complex<float>>::const_contiter inPtr =
      data.cbegin() + startChannel * polCount;
  const size_t selectedChannelCount = endChannel - startChannel;

  if (polOut == aocommon::Polarization::Instrumental) {
    if (polsIn.size() != 4) {
      throw std::runtime_error(
          "This mode requires the four polarizations to be present in the "
          "measurement set");
    }
    for (size_t ch = 0; ch != selectedChannelCount * polsIn.size(); ++ch) {
      if (IsCFinite(*inPtr))
        dest[ch] = *inPtr;
      else
        dest[ch] = 0;
      ++inPtr;
    }
  } else if (polOut == aocommon::Polarization::DiagonalInstrumental) {
    if (polsIn.size() == 4) {
      size_t ch = 0;
      while (ch != selectedChannelCount * 2) {
        if (IsCFinite(*inPtr))
          dest[ch] = *inPtr;
        else
          dest[ch] = 0;
        inPtr += 3;  // jump from xx to yy
        ++ch;
        if (IsCFinite(*inPtr))
          dest[ch] = *inPtr;
        else
          dest[ch] = 0;
        ++inPtr;
        ++ch;
      }
    } else if (polsIn.size() == 2) {
      for (size_t ch = 0; ch != selectedChannelCount * 2; ++ch) {
        if (IsCFinite(*inPtr))
          dest[ch] = *inPtr;
        else
          dest[ch] = 0;
        ++inPtr;
      }
    } else
      throw std::runtime_error(
          "Diagonal instrument visibilities requested, but this requires 2 or "
          "4 polarizations in the data");
  } else if (size_t polIndex;
             aocommon::Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
    inPtr += polIndex;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      if (IsCFinite(*inPtr))
        dest[ch] = *inPtr;
      else
        dest[ch] = 0;
      inPtr += polCount;
    }
  } else {
    // Copy the right visibilities with conversion if necessary.
    switch (polOut) {
      case aocommon::Polarization::StokesI: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (!hasXX || !hasYY) {
          bool hasRR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, polsIn, polIndexA);
          bool hasLL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, polsIn, polIndexB);
          if (!hasRR || !hasLL)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes I) from available "
                "polarizations");
        }

        for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
          inPtr += polIndexA;
          casacore::Complex val = *inPtr;
          inPtr += polIndexB - polIndexA;

          // I = (XX + YY) / 2
          val = (*inPtr + val) * 0.5f;

          if (IsCFinite(val))
            dest[ch] = val;
          else
            dest[ch] = 0.0;

          inPtr += polCount - polIndexB;
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (hasXX && hasYY) {
          // Convert to StokesQ from XX and YY
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // Q = (XX - YY)/2
            val = (val - *inPtr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        } else {
          // Convert to StokesQ from RR and LL
          bool hasRL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, polsIn, polIndexA);
          bool hasLR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, polsIn, polIndexB);
          if (!hasRL || !hasLR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes Q) from available "
                "polarizations");
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // Q = (RL + LR)/2
            val = (*inPtr + val) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        }
      } break;
      case aocommon::Polarization::StokesU: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (hasXY && hasYX) {
          // Convert to StokesU from XY and YX
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // U = (XY + YX)/2
            val = (val + *inPtr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        } else {
          // Convert to StokesU from RR and LL
          bool hasRL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RL, polsIn, polIndexA);
          bool hasLR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LR, polsIn, polIndexB);
          if (!hasRL || !hasLR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes U) from available "
                "polarizations");
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // U = -i (RL - LR)/2
            val = (val - *inPtr) * 0.5f;
            val = casacore::Complex(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        }
      } break;
      case aocommon::Polarization::StokesV: {
        size_t polIndexA = 0, polIndexB = 0;
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (hasXY && hasYX) {
          // Convert to StokesV from XX and YY
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // V = -i(XY - YX)/2
            val = (val - *inPtr) * 0.5f;
            val = casacore::Complex(val.imag(), -val.real());

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        } else {
          // Convert to StokesV from RR and LL
          bool hasRL = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::RR, polsIn, polIndexA);
          bool hasLR = aocommon::Polarization::TypeToIndex(
              aocommon::Polarization::LL, polsIn, polIndexB);
          if (!hasRL || !hasLR)
            throw std::runtime_error(
                "Can not form requested polarization (Stokes V) from available "
                "polarizations");
          for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
            inPtr += polIndexA;
            casacore::Complex val = *inPtr;
            inPtr += polIndexB - polIndexA;

            // V = (RR - LL)/2
            val = (val - *inPtr) * 0.5f;

            if (IsCFinite(val))
              dest[ch] = val;
            else
              dest[ch] = 0.0;

            inPtr += polCount - polIndexB;
          }
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
    }
  }
}

template <typename NumType>
void PartitionClass::CopyWeights(
    NumType* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut) {
  const size_t polCount = polsIn.size();
  casacore::Array<std::complex<float>>::const_contiter dataPtr =
      data.cbegin() + startChannel * polCount;
  casacore::Array<float>::const_contiter weightPtr =
      weights.cbegin() + startChannel * polCount;
  casacore::Array<bool>::const_contiter flagPtr =
      flags.cbegin() + startChannel * polCount;
  const size_t selectedChannelCount = endChannel - startChannel;

  size_t polIndex;
  if (polOut == aocommon::Polarization::Instrumental) {
    for (size_t ch = 0; ch != selectedChannelCount * polsIn.size(); ++ch) {
      if (!*flagPtr && IsCFinite(*dataPtr))
        // The factor of 4 is to be consistent with StokesI
        // It is for having conjugate visibilities and because IDG doesn't
        // separately count XX and YY visibilities
        dest[ch] = *weightPtr * 4.0f;
      else
        dest[ch] = 0.0f;
      dataPtr++;
      weightPtr++;
      flagPtr++;
    }
  } else if (polOut == aocommon::Polarization::DiagonalInstrumental) {
    if (polsIn.size() == 4) {
      size_t ch = 0;
      while (ch != selectedChannelCount * 2) {
        if (!*flagPtr && IsCFinite(*dataPtr))
          // See explanation above for factor of 4
          dest[ch] = *weightPtr * 4.0f;
        else
          dest[ch] = 0.0f;
        dataPtr += 3;  // jump from xx to yy
        weightPtr += 3;
        flagPtr += 3;
        ++ch;
        if (!*flagPtr && IsCFinite(*dataPtr))
          dest[ch] = *weightPtr * 4.0f;
        else
          dest[ch] = 0.0f;
        ++dataPtr;
        ++weightPtr;
        ++flagPtr;
        ++ch;
      }
    } else if (polsIn.size() == 2) {
      for (size_t ch = 0; ch != selectedChannelCount * 2; ++ch) {
        if (!*flagPtr && IsCFinite(*dataPtr))
          dest[ch] = *weightPtr * 4.0f;
        else
          dest[ch] = 0.0f;
        ++dataPtr;
        ++weightPtr;
        ++flagPtr;
      }
    }
  } else if (aocommon::Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
    dataPtr += polIndex;
    weightPtr += polIndex;
    flagPtr += polIndex;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      if (!*flagPtr && IsCFinite(*dataPtr))
        dest[ch] = *weightPtr;
      else
        dest[ch] = 0.0f;
      dataPtr += polCount;
      weightPtr += polCount;
      flagPtr += polCount;
    }
  } else {
    size_t polIndexA = 0, polIndexB = 0;
    switch (polOut) {
      case aocommon::Polarization::StokesI: {
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (!hasXY || !hasYX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsIn, polIndexB);
        }
      } break;
      case aocommon::Polarization::StokesQ: {
        bool hasXX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XX, polsIn, polIndexA);
        bool hasYY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YY, polsIn, polIndexB);
        if (!hasXX || !hasYY) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              polsIn, polIndexB);
        }
      } break;
      case aocommon::Polarization::StokesU: {
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (!hasXY || !hasYX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RL,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LR,
                                              polsIn, polIndexB);
        }
      } break;
      case aocommon::Polarization::StokesV: {
        bool hasXY = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::XY, polsIn, polIndexA);
        bool hasYX = aocommon::Polarization::TypeToIndex(
            aocommon::Polarization::YX, polsIn, polIndexB);
        if (!hasXY || !hasYX) {
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::RR,
                                              polsIn, polIndexA);
          aocommon::Polarization::TypeToIndex(aocommon::Polarization::LL,
                                              polsIn, polIndexB);
        }
      } break;
      default:
        throw std::runtime_error(
            "Could not convert ms polarizations to requested polarization");
        break;
    }

    weightPtr += polIndexA;
    dataPtr += polIndexA;
    flagPtr += polIndexA;
    for (size_t ch = 0; ch != selectedChannelCount; ++ch) {
      NumType w;
      if (!*flagPtr && IsCFinite(*dataPtr))
        w = *weightPtr * 4.0f;
      else
        w = 0.0f;
      dataPtr += polIndexB - polIndexA;
      weightPtr += polIndexB - polIndexA;
      flagPtr += polIndexB - polIndexA;
      if (!*flagPtr && IsCFinite(*dataPtr))
        w = std::min<NumType>(w, *weightPtr * 4.0f);
      else
        w = 0.0f;
      dest[ch] = w;
      weightPtr += polCount - polIndexB + polIndexA;
      dataPtr += polCount - polIndexB + polIndexA;
      flagPtr += polCount - polIndexB + polIndexA;
    }
  }
}

template void PartitionClass::CopyWeights<float>(
    float* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut);

template void PartitionClass::CopyWeights<std::complex<float>>(
    std::complex<float>* dest, size_t startChannel, size_t endChannel,
    const std::vector<aocommon::PolarizationEnum>& polsIn,
    const casacore::Array<std::complex<float>>& data,
    const casacore::Array<float>& weights, const casacore::Array<bool>& flags,
    aocommon::PolarizationEnum polOut);

void PartitionClass::preprocessPartition(const string& dataColumnName,
                                         const Settings& settings) {
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
  polarizationsPerFile =
      aocommon::Polarization::GetVisibilityCount(*polsOut.begin());
  temporaryDirectory = settings.temporaryDirectory;

  channelParts = channels.size();
  if (channelParts != 1) {
    Logger::Debug << "Partitioning in " << channels.size() << " channels:";
    for (size_t i = 0; i != channels.size(); ++i)
      Logger::Debug << ' ' << channels[i].dataDescId << ':' << channels[i].start
                    << '-' << channels[i].end;
  }
  Logger::Debug << '\n';

  // Ordered as files[pol x channelpart]
  files.resize(channelParts * polsOut.size());
  max_channels = GetMaxChannels(channels);

  // Each data desc id needs a separate meta file because they can have
  // different uvws and other info.
  size_t fileIndex = 0;
  for (size_t part = 0; part != channelParts; ++part) {
    for (aocommon::PolarizationEnum p : polsOut) {
      PartitionFiles& f = files[fileIndex];
      std::string partPrefix = getPartPrefix(
          msOutPath, part, p, channels[part].dataDescId, temporaryDirectory);
      f.data = std::make_unique<std::ofstream>(partPrefix + ".tmp");
      f.weight = std::make_unique<std::ofstream>(partPrefix + "-w.tmp");
      if (initialModelRequired)
        f.model = std::make_unique<std::ofstream>(partPrefix + "-m.tmp");
      f.data->seekp(PartHeader::BINARY_SIZE, std::ios::beg);

      ++fileIndex;
    }
  }

  // This maps dataDescId to spw index.
  selectedDataDescIds = getDataDescIdMap(channels);
  if (settings.parallelReordering == 1)
    Logger::Info << "Reordering " << msPath << " into " << channelParts << " x "
                 << polsOut.size() << " parts.\n";

  // Write header of meta file, one meta file for each data desc id
  // TODO rather than writing we can just skip and write later
  metaFiles.resize(selectedDataDescIds.size());
  for (const std::pair<const size_t, size_t>& p : selectedDataDescIds) {
    const size_t dataDescId = p.first;
    const size_t spwIndex = p.second;
    std::string metaFilename =
        getMetaFilename(msOutPath, temporaryDirectory, dataDescId);
    metaFiles[spwIndex] = std::make_unique<std::ofstream>(metaFilename);
    MetaHeader metaHeader;
    metaHeader.selectedRowCount = 0;  // not yet known
    metaHeader.filenameLength = msOutPath.size();
    metaHeader.startTime = 0;  // rowProvider->StartTime();
    metaHeader.Write(*metaFiles[spwIndex]);
    metaFiles[spwIndex]->write(msOutPath.c_str(), msOutPath.size());
    if (!metaFiles[spwIndex]->good())
      throw std::runtime_error("Error writing to temporary file " +
                               metaFilename);
  }
}

void PartitionClass::processPartition(dp3::base::DPBuffer* buffer,
                                      const Settings& settings) {
  casacore::MeasurementSet msObj(msPath);
  const std::map<size_t, std::vector<aocommon::PolarizationEnum>>
      msPolarizationsPerDataDescId =
          GetMSPolarizationsPerDataDescId(channels, msObj);

  // Check for simulated noise.
  std::normal_distribution<float> _distribution(0.0, 1.0);
  if (settings.baselineDependentAveragingInWavelengths == 0.0) {
    if (settings.simulateNoise) {
      if (settings.simulatedBaselineNoiseFilename.empty())
        _distribution =
            std::normal_distribution<float>(0.0, settings.simulatedNoiseStdDev);
      else {
        std::ifstream file(settings.simulatedBaselineNoiseFilename);
        if (!file)
          throw std::runtime_error("Can't open baseline noise file " +
                                   settings.simulatedBaselineNoiseFilename);
        _noiseMap = NoiseMap(file);
      }
    }
  }

  // Write actual data
  dataBuffer.assign(max_channels * polarizationsPerFile, 0.0);
  std::vector<float> weightBuffer(polarizationsPerFile * max_channels);

  // Get DP3 time frame details.
  unsigned int n_baselines = buffer->GetFlags().shape(0);
  unsigned int n_channels = buffer->GetFlags().shape(1);
  unsigned int n_correlations = buffer->GetFlags().shape(2);

  dp3::base::DPBuffer::FlagsType buffFlags = buffer->GetFlags();
  dp3::base::DPBuffer::UvwType buffUvw = buffer->GetUvw();
  dp3::base::DPBuffer::WeightsType buffWeights = buffer->GetWeights();
  dp3::base::DPBuffer::DataType buffData = buffer->GetData();
  casacore::Vector<dp3::common::rownr_t> rowNum = buffer->GetRowNumbers();
  std::vector<int> antenna1List = infoObj.getAnt1();
  std::vector<int> antenna2List = infoObj.getAnt2();

  casacore::ScalarColumn<casacore::Int> dataDescTable(
      msObj, casacore::MS::columnName(casacore::MS::DATA_DESC_ID));
  casacore::Vector<long long unsigned int> rowNumArray =
      buffer->GetRowNumbers();
  // uint32_t dataDescId = infoObj.spectralWindow();

  // Casacore arrays used to write the meta, data, weight and model files.
  casacore::Array<std::complex<float>> dataArray;
  casacore::Array<std::complex<float>> modelArray;
  casacore::Array<float> weightSpectrumArray;
  casacore::Array<bool> flagArray;

  static size_t selectedRowsTotal = 0;
  selectedRowCountPerSpwIndex.resize(selectedDataDescIds.size(), 0);

  for (size_t bl = 0; bl < n_baselines; bl++) {
    MetaRecord meta;
    size_t strideVal = n_channels * n_correlations;
    bool* flagPtr = &buffFlags(bl, 0, 0);
    float* weightPtr = &buffWeights(bl, 0, 0);
    std::complex<float>* dataPtr = &buffData(bl, 0, 0);

    flagArray = casacore::Vector<bool>(flagPtr, strideVal, 0);
    weightSpectrumArray = casacore::Vector<float>(weightPtr, strideVal, 0);

    if (settings.simulateNoise) {
      std::mt19937 _rng;
      const float stddev =
          _noiseMap.Empty()
              ? 1.0
              : _noiseMap.GetNoiseValue(antenna1List[bl], antenna2List[bl]);

      casacore::Array<std::complex<float>>::contiter iter;
      for (iter = dataArray.cbegin(); iter != dataArray.cend(); ++iter) {
        if (std::isfinite(iter->real()) && std::isfinite(iter->imag())) {
          iter->real(_distribution(_rng) * stddev);
          iter->imag(_distribution(_rng) * stddev);
        } else {
          iter->real(std::numeric_limits<float>::quiet_NaN());
          iter->imag(std::numeric_limits<float>::quiet_NaN());
        }
      }
    } else {
      dataArray = casacore::Vector<std::complex<float>>(dataPtr, strideVal, 0);
    }

    meta.u = buffUvw(bl, 0);
    meta.v = buffUvw(bl, 1);
    meta.w = buffUvw(bl, 2);
    meta.antenna1 = antenna1List[bl];
    meta.antenna2 = antenna2List[bl];
    meta.fieldId = 0;
    meta.time = buffer->GetTime();

    uint64_t dataDescId = dataDescTable.get(rowNumArray[bl]);
    const size_t spwIndex = selectedDataDescIds.find(dataDescId)->second;
    ++selectedRowCountPerSpwIndex[spwIndex];
    ++selectedRowsTotal;
    std::ofstream& metaFile = *metaFiles[spwIndex];
    meta.Write(metaFile);
    if (!metaFile.good())
      throw std::runtime_error("Error writing to temporary file");

    if (initialModelRequired)
      modelArray = casacore::Vector<std::complex<float>>(dataPtr, strideVal, 0);

    size_t fileIndex = 0;
    for (size_t part = 0; part != channelParts; ++part) {
      if (channels[part].dataDescId == int(dataDescId)) {
        const size_t partStartCh = channels[part].start;
        const size_t partEndCh = channels[part].end;
        const std::vector<aocommon::PolarizationEnum>& msPolarizations =
            msPolarizationsPerDataDescId.find(dataDescId)->second;

        for (aocommon::PolarizationEnum p : polsOut) {
          PartitionFiles& f = files[fileIndex];
          CopyData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations,
                   dataArray, p);
          f.data->write(reinterpret_cast<char*>(dataBuffer.data()),
                        (partEndCh - partStartCh) *
                            sizeof(std::complex<float>) * polarizationsPerFile);
          if (!f.data->good())
            throw std::runtime_error("Error writing to temporary data file");

          if (initialModelRequired) {
            CopyData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations,
                     modelArray, p);
            f.model->write(reinterpret_cast<char*>(dataBuffer.data()),
                           (partEndCh - partStartCh) *
                               sizeof(std::complex<float>) *
                               polarizationsPerFile);
            if (!f.model->good())
              throw std::runtime_error(
                  "Error writing to temporary model data file");
          }

          CopyWeights(weightBuffer.data(), partStartCh, partEndCh,
                      msPolarizations, dataArray, weightSpectrumArray,
                      flagArray, p);
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
  }
  // Logger::Info << "Total selected rows: " << selectedRowsTotal << '\n';
}

void PartitionClass::postprocess() {
  // Rewrite meta headers to include selected row count
  for (const std::pair<const size_t, size_t>& p : selectedDataDescIds) {
    const size_t spwIndex = p.second;
    MetaHeader metaHeader;
    metaHeader.selectedRowCount = selectedRowCountPerSpwIndex[spwIndex];
    metaHeader.filenameLength = msOutPath.size();
    metaHeader.startTime = 0;  // rowProvider->StartTime();
    metaFiles[spwIndex]->seekp(0);
    metaHeader.Write(*metaFiles[spwIndex]);
    metaFiles[spwIndex]->write(msOutPath.c_str(), msOutPath.size());
  }

  // Write header to parts and write empty model files (if requested)
  PartHeader header;
  header.hasModel = includeModel;
  size_t fileIndex = 0;
  dataBuffer.assign(max_channels * polarizationsPerFile, 0.0);

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
            msOutPath, part, *p, header.dataDescId, temporaryDirectory);
        std::ofstream modelFile(partPrefix + "-m.tmp");
        const size_t selectedRowCount = selectedRowCountPerSpwIndex
            [selectedDataDescIds.find(channels[part].dataDescId)->second];
        for (size_t i = 0; i != selectedRowCount; ++i) {
          modelFile.write(reinterpret_cast<char*>(dataBuffer.data()),
                          header.channelCount * sizeof(std::complex<float>) *
                              polarizationsPerFile);
        }
      }
    }
  }
}