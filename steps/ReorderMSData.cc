
#include "ReorderMSData.h"

#include <limits>


using dp3::common::ParameterSet;
using dp3::base::DPInfo;

namespace dp3 {
    namespace steps {

        /// Construct the object.
        /// Parameters are obtained from the parset using the given prefix.
        Reorder::Reorder(const common::ParameterSet& parset, const string& prefix):settingsObj(cleanObj.GetSettings())
        {
            kRequiredFields = kDataField | kFlagsField | kWeightsField | kUvwField;
            kProvidedFields = kRequiredFields;
            configureSettings();
            
            cleanObj.RunClean();
            // cleanObj.RunPredict();
            cleanObj.performReordering(predictMode);
        }

        void Reorder::configureSettings()
        {
            predictMode = false; // parset.getString("predictMode");
            cleanObj._settings.filenames.push_back(parset.getString("msin"));
            // Other Settings
        }

        /// Process the data.
        /// It keeps the data.
        /// When processed, it invokes the process function of the next step.
        bool Reorder::process(std::unique_ptr<base::DPBuffer> buffer)
        {
            cleanObj.partitionObj.processPartition(buffer.get());
            getNextStep()->process(std::move(buffer));
            return true;
        }

        /// Show the step parameters.
        void Reorder::show(std::ostream& outputStream) const
        {
            outputStream << "Test output" << std::endl;
        } 

        void Reorder::finish() { 
            cleanObj.partitionObj.postprocess();
            getNextStep()->finish(); 
        }

    }  // namespace steps
}