
#ifndef DP3_STEPS_REORDER_H_
#define DP3_STEPS_REORDER_H_

#include <memory>

#include <aocommon/staticfor.h>

#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include <dp3/steps/Step.h>
#include "../common/Timer.h"
#include "../common/ParameterSet.h"
#include <iostream>
#include <string>

#include "reorderHelperFiles/include/wsclean.h"

namespace dp3 {

namespace common {
class ParameterSet;
class ParameterValue;
}  // namespace common

namespace steps {

class Reorder : public Step {
 public:
  /// Required and provided fields are fixed. These constants simplify
  /// get*Fields implementations for steps that have Reorder as sub-step.
  common::Fields kRequiredFields;
  common::Fields kProvidedFields;

  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  Reorder(const common::ParameterSet& parset, const string& prefix);
  
  ~Reorder() override{}

  common::Fields getRequiredFields() const override { return kRequiredFields; }

  common::Fields getProvidedFields() const override { return kProvidedFields; }

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  bool process(std::unique_ptr<base::BDABuffer> buffer) override {
    // Do the BDA processing here.
    return getNextStep()->process(std::move(buffer));
  }
  
  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo& infoObjIn) override { 
    std::cout << "Updated Info Object" << std::endl;
    cleanObj.partitionObj.infoObj = infoObjIn;
    Step::updateInfo(infoObjIn); 
  }

  /// Show the step parameters.
  void show(std::ostream& outputStream) const override;
  
  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override {}

  /// Get the value in Hertz of a string like "1000 MHz". If unit is
  /// omitted it defaults to Hertz
  static double getFreqHz(const string& freqstr) { return 0; }

 private:
  /// Update itsBuf so it contains averages.
  bool predictMode;
  WSClean cleanObj;
  Settings& settingsObj;
  void getSettings(const common::ParameterSet& parset);
    
};

}  // namespace steps
}  // namespace dp3

#endif
