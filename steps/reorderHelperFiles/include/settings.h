#ifndef WSCLEAN_SETTINGS_H
#define WSCLEAN_SETTINGS_H

#include <cassert>
#include <optional>

#include "weightmode.h"
#include "msselection.h"

#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>

enum class DirectFTPrecision { Float, Double, LongDouble };
enum class GridderType { WStacking, WGridder, TunedWGridder, DirectFT, IDG };

/**
 * This class describes all settings for a single WSClean run.
 * @sa WSClean
 */
class Settings {
 public:
  Settings();

  void Validate() const;
  void checkPolarizations() const;

  size_t GetFeatherSize() const {
    if (featherSize) {
      return *featherSize;
    } else {
      // Return the default: 1% of sqrt(width * height)
      return std::ceil(std::sqrt(trimmedImageWidth * trimmedImageHeight) *
                       0.01);
    }
  }

  std::vector<std::string> filenames;
  GridderType gridderType = GridderType::WGridder;
  bool continuedRun;
  size_t parallelReordering;

  /// Maps output channel indices to node indices, when using MPI.
  std::vector<size_t> fieldIds;
  std::set<aocommon::PolarizationEnum> polarizations;
  std::set<size_t> spectralWindows;
  std::string temporaryDirectory;
  std::string simulatedBaselineNoiseFilename;
  std::string atermConfigFilename;
  size_t ddPsfGridHeight, ddPsfGridWidth;
  bool subtractModel, modelUpdateRequired;
  bool diagonalSolutions = false;
  bool gridWithBeam;
  bool simulateNoise;
  double baselineDependentAveragingInWavelengths;
  double simulatedNoiseStdDev;
  double deconvolutionMGain;

  size_t startTimestep;
  size_t endTimestep;
  double minUVWInMeters;
  double maxUVWInMeters;
  enum MSSelection::EvenOddSelection evenOddTimesteps;

  bool hasShift;
  bool joinedFrequencyDeconvolution;
  bool joinedPolarizationDeconvolution;
  bool divideChannelsByGaps;
  aocommon::UVector<double> divideChannelFrequencies;
  std::optional<size_t> featherSize;
  std::vector<float> spectralCorrection;
  std::string facetRegionFilename;
  std::string dataColumnName;
  size_t trimmedImageWidth, trimmedImageHeight;
  size_t deconvolutionIterationCount;
  size_t intervalsOut;
  size_t endChannel, startChannel;
  size_t channelsOut;
  double pixelScaleX;
  double pixelScaleY;
  double imagePadding;
  double shiftRA;
  double shiftDec;
  double spectralCorrectionFrequency;

  MSSelection GetMSSelection() const {
    MSSelection selection;
    selection.SetInterval(startTimestep, endTimestep);
    selection.SetFieldIds(fieldIds);
    selection.SetMinUVWInM(minUVWInMeters);
    selection.SetMaxUVWInM(maxUVWInMeters);
    selection.SetEvenOrOddTimesteps(evenOddTimesteps);
    return selection;
  }
};

inline Settings::Settings()
    : filenames(),
      gridderType(GridderType::WGridder),
      continuedRun(false),
      parallelReordering(4),

      fieldIds({0}),
      polarizations({aocommon::Polarization::StokesI}),
      temporaryDirectory(),
      simulatedBaselineNoiseFilename(),
      atermConfigFilename(),
      ddPsfGridHeight(1),
      ddPsfGridWidth(1),
      subtractModel(false),
      modelUpdateRequired(true),
      diagonalSolutions(false),
      gridWithBeam(false),
      simulateNoise(false),
      baselineDependentAveragingInWavelengths(0.0),
      simulatedNoiseStdDev(0.0),
      deconvolutionMGain(1.0),

      hasShift(false),
      joinedFrequencyDeconvolution(false),
      joinedPolarizationDeconvolution(false),
      divideChannelsByGaps(false),
      divideChannelFrequencies(),
      spectralCorrection(),
      facetRegionFilename(),
      dataColumnName(),
      trimmedImageWidth(0),
      trimmedImageHeight(0),
      deconvolutionIterationCount(0),
      intervalsOut(1),
      endChannel(0),
      startChannel(0),
      channelsOut(1),
      pixelScaleX(0.0),
      pixelScaleY(0.0),
      imagePadding(1.2),
      shiftRA(0.0),
      shiftDec(0.0),
      spectralCorrectionFrequency(0.0) {}

#endif
