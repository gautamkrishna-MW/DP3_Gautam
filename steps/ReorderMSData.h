
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

  /// Process the data.
  /// It keeps the data.
  /// When processed, it invokes the process function of the next step.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  common::Fields getRequiredFields() const override { 
    return kRequiredFields; 
  }

  common::Fields getProvidedFields() const override { 
    return kProvidedFields; 
  }

  bool process(std::unique_ptr<base::BDABuffer> buffer) override {
    // Do the BDA processing here.
    return getNextStep()->process(std::move(buffer));
  }
  
  /// Update the general info.
  void updateInfo(const base::DPInfo& infoObjIn) override { 
    // std::cout << "Updated Info Object" << std::endl;
    cleanObj.partitionObj.infoObj = infoObjIn;
    Step::updateInfo(infoObjIn); 
  }

  void configureSettings(const common::ParameterSet& parset);
  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override {}
  /// Show the step parameters.
  void show(std::ostream& outputStream) const override;
  /// Finish the processing of this step and subsequent steps.
  void finish() override;

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
