
#include "ReorderMSData.h"

#include <limits>


using dp3::common::ParameterSet;
using dp3::base::DPInfo;

namespace dp3 {
    namespace steps {

        void Reorder::getSettings(const common::ParameterSet& parset)
        {
            settingsObj.filenames.push_back(parset.getString("msin"));

            // isPredictMode
            settingsObj.deconvolutionMGain = 1.0;
            settingsObj.subtractModel = false;
            settingsObj.continuedRun = false;
            settingsObj.parallelReordering = 4;
            settingsObj.modelUpdateRequired = true;
            settingsObj.gridderType = GridderType::WGridder;
            settingsObj.polarizations = {aocommon::Polarization::StokesI};
            settingsObj.ddPsfGridWidth = 1;
            settingsObj.ddPsfGridHeight = 1;
            settingsObj.gridWithBeam = false;
            settingsObj.diagonalSolutions = false;
            settingsObj.temporaryDirectory= std::string("");
            settingsObj.fieldIds = {0};
            settingsObj.baselineDependentAveragingInWavelengths = 0.0;
            settingsObj.simulatedBaselineNoiseFilename = std::string("");
            settingsObj.simulatedNoiseStdDev = 0.0;
            settingsObj.simulateNoise = false;
        }

        /// Construct the object.
        /// Parameters are obtained from the parset using the given prefix.
        Reorder::Reorder(const common::ParameterSet& parset, const string& prefix):settingsObj(cleanObj.GetSettings())
        {
            kRequiredFields = kDataField | kFlagsField | kWeightsField | kUvwField;
            kProvidedFields = kRequiredFields;
            std::cout << "Reorder Call: " << prefix << std::endl;
            std::cout << "MS Path: " << parset.getString("msin") << std::endl;
            std::cout << "MS Path: " << parset.getString("msout") << std::endl;
            settingsObj = cleanObj.GetSettings();

            // Function to get the settings of reordering.
            getSettings(parset);
            predictMode = false;
            cleanObj.RunClean();
            // cleanObj.RunPredict();
            cleanObj.performReordering(predictMode);
        }

        /// Process the data.
        /// It keeps the data.
        /// When processed, it invokes the process function of the next step.
        bool Reorder::process(std::unique_ptr<base::DPBuffer> buffer)
        {
            // static int counter = 0;
            // dp3::base::DPBuffer::FlagsType buffFlags = buffer->GetFlags();

            // static int newCount = 0;
            // std::cout << "CounterTop: " << newCount++ << std::endl;
            // std::cout << "Time: " << buffer->GetTime() << std::endl;

            // for (auto& row : buffer->GetRowNumbers())
            // {
            //     std::cout << "\n\nRow: " << row << std::endl;
            //     std::cout << "Ant: " << getInfo().getAnt1()[row] << ", " << getInfo().getAnt2()[row] << std::endl;
            //     std::cout << "SPW: " << getInfo().spectralWindow() << std::endl;
            //     std::cout << "Corr: " << getInfo().ncorr() << std::endl;
            //     std::cout << "Chan: " << getInfo().nchan() << std::endl;
            //     std::cout << "BasL: " << getInfo().nbaselines() << std::endl;
            //     std::cout << "Time: " << getInfo().ntime() << std::endl;
            //     std::cout << "nAnt: " << getInfo().nantenna() << std::endl;
            // }

            // for (auto rows : buffer->GetRowNumbers())
            //     std::cout << "Sol: " << buffer->GetFlags()[] << std::endl;
            // std::cout << "Buff size: " << buffer->GetFlags().size() << std::endl;
            // std::unique_ptr<base::DPBuffer> newBuffer = std::move(buffer);
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