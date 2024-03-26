#include "../include/settings.h"

#include "../include/facetreader.h"

#include <aocommon/logger.h>

#include <schaapcommon/h5parm/h5parm.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <sstream>

using aocommon::Logger;

void Settings::Validate() const {
  if (diagonalSolutions) {
    if (polarizations.size() != 1 ||
        *polarizations.begin() != aocommon::PolarizationEnum::StokesI) {
      throw std::runtime_error(
          "-diagonal-solutions can only be used when making Stokes I images");
    }
  }

  if (gridderType == GridderType::IDG) {
    const bool stokesIOnly =
        polarizations.size() == 1 &&
        *polarizations.begin() == aocommon::Polarization::StokesI;
    const bool allStokes =
        aocommon::Polarization::HasFullStokesPolarization(polarizations) &&
        polarizations.size() == 4;

    if (!allStokes && !stokesIOnly) {
      throw std::runtime_error(
          "When using IDG, it is only possible to either image Stokes I or to "
          "image all 4 Stokes polarizations: use -pol i or -pol iquv.");
    }

    if (gridderType != GridderType::IDG && !atermConfigFilename.empty())
      throw std::runtime_error(
          "Use of an aterm config file required IDG enabled: add -use-idg");

    if (baselineDependentAveragingInWavelengths != 0.0) {
      throw std::runtime_error(
          "IDG cannot be combined with (internally computed) "
          "baseline-dependent averaging. Please remove baseline-averaging "
          "option from your command.");
    }

    for (const auto& filename : filenames) {
      casacore::MeasurementSet ms(filename);
      const std::string& bda_factors = "BDA_FACTORS";
      const bool has_bda = ms.keywordSet().isDefined(bda_factors) &&
                           (ms.keywordSet().asTable(bda_factors).nrow() > 0);
      if (has_bda) {
        throw std::runtime_error(
            "IDG cannot be combined with the baseline-dependently averaged "
            "measurement set " +
            filename);
      }
    }
  }

  if (gridderType != GridderType::IDG && !atermConfigFilename.empty())
    throw std::runtime_error(
        "Use of an aterm config file required IDG enabled: add -use-idg");

  if (gridWithBeam && !atermConfigFilename.empty())
    throw std::runtime_error(
        "Use of an aterm config file can't be combined with -grid-with-beam: "
        "add the beam to your aterm config and remove -grid-with-beam from the "
        "command line");

  if (gridWithBeam && gridderType != GridderType::IDG)
    throw std::runtime_error(
        "Can't grid with the beam without IDG: specify '-use-idg' to use IDG.");

  if (baselineDependentAveragingInWavelengths != 0.0) {
    if (modelUpdateRequired)
      throw std::runtime_error(
          "Baseline dependent averaging can not update the model column (yet) "
          "-- you have to add -no-update-model-required.");
  }

  checkPolarizations();
}

void Settings::checkPolarizations() const {
  bool hasXY = polarizations.count(aocommon::Polarization::XY) != 0;
  bool hasYX = polarizations.count(aocommon::Polarization::YX) != 0;

  if ((hasXY && !hasYX) || (!hasXY && hasYX))
    throw std::runtime_error(
        "You are imaging only one of the XY or YX polarizations. This is not "
        "possible -- you have to specify both XY and YX polarizations (the "
        "output of imaging both polarizations will be the XY and imaginary XY "
        "images).");
}
