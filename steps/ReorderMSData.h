
#ifndef DP3_STEPS_REORDER_H_
#define DP3_STEPS_REORDER_H_

#include <memory>

#include <aocommon/staticfor.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/steps/Step.h>
#include "../common/Timer.h"
#include <iostream>

#include "reorderHelperFiles/include/wsclean.h"

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

class Reorder : public Step {
 public:
  /// Required and provided fields are fixed. These constants simplify
  /// get*Fields implementations for steps that have Averagers as sub-step.
  common::Fields kRequiredFields;
  common::Fields kProvidedFields;

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Reorder(const common::ParameterSet& parset, const string& prefix)
  {
    kRequiredFields = kDataField | kFlagsField | kWeightsField | kUvwField;
    kProvidedFields = kRequiredFields;
    std::cout << "Reorder Call: " << prefix << std::endl;
    std::cout << "MS Path: " << parset.getString("msin") << std::endl;
    std::cout << "MS Path: " << parset.getString("msout") << std::endl;

    WSClean cleanObj;
    Settings& settingsObj = cleanObj.GetSettings();
    settingsObj.filenames.push_back(parset.getString("msin"));
    // cleanObj.RunClean();
    cleanObj.RunPredict();

  }

  /// Construct the object using the given parameters.
  Reorder(const string& stepname, unsigned int nchanAvg,
           unsigned int ntimeAvg){}
  Reorder(const string& stepname, double freq_resolution,
           double time_resolution){}

  ~Reorder() override{}

  common::Fields getRequiredFields() const override { return kRequiredFields; }

  common::Fields getProvidedFields() const override { return kProvidedFields; }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override
  {}

  /// Finish the processing of this step and subsequent steps.
  void finish() override
  {}

  /// Update the general info.
  void updateInfo(const base::DPInfo& infoObj) override
  {}

  /// Show the step parameters.
  void show(std::ostream& outputStream) const override
  {
    outputStream << "Test output" << std::endl;
  }

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override
  {}

  /// Get the value in Hertz of a string like "1000 MHz". If unit is
  /// omitted it defaults to Hertz
  static double getFreqHz(const string& freqstr)
  {}

 private:
  /// Update itsBuf so it contains averages.
  
};

}  // namespace steps
}  // namespace dp3

#endif
