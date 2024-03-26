
#include "ReorderMSData.h"

#include <limits>
#include <iomanip>


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

            // Configure Settings
            std::string inputMSPath(parset.getString("msin"));
            std::string outputMSPath(parset.getString("msout"));

            cleanObj.partitionObj.msOutPath = outputMSPath;
            if (outputMSPath.empty())
                cleanObj.partitionObj.msOutPath = inputMSPath;
            
            cleanObj._settings.filenames.push_back(inputMSPath);
            predictMode = parset.getBool("predictMode", false);
            configureSettings(parset);
            
            if (!predictMode)
                cleanObj.RunClean();
            else
                cleanObj.RunPredict();

            cleanObj.performReordering(predictMode);
        }

        /// Process the data.
        /// It keeps the data.
        /// When processed, it invokes the process function of the next step.
        bool Reorder::process(std::unique_ptr<base::DPBuffer> buffer)
        {
            // static size_t numTimes = 0;
            // std::cout << "Time Slot: " << numTimes++ << std::endl;
            // std::cout << "Time: " << std::setprecision(10) << buffer->GetTime() << std::endl;
            // std::cout << "Num BL: " << buffer->GetFlags().shape(0) << std::endl;
            // std::cout << "Num Rows: " << buffer->GetRowNumbers().size() << std::endl;
            // std::cout << "Num UVW: " << buffer->GetUvw().shape(0) << ", " << buffer->GetUvw().shape(1) << std::endl;
            // std::cout << "Num Data: " << buffer->GetData().shape(0) << ", " << buffer->GetData().shape(1) << ", " << buffer->GetData().shape(2) << std::endl;
            // std::cout << "Num Weight: " << buffer->GetWeights().shape(0) << ", " << buffer->GetWeights().shape(1) << ", " << buffer->GetWeights().shape(2) << std::endl;
            // std::cout << "Process done" << std::endl;
            // std::cout << "\n\n" << std::endl;
            cleanObj.partitionObj.processPartition(buffer.get(), cleanObj._settings);
            getNextStep()->process(std::move(buffer));
            return true;
        }

        /// Show the step parameters.
        void Reorder::show(std::ostream& os) const
        {
            os << "Reorder\n";
            os << "  input MS:       " << cleanObj.partitionObj.msPath << '\n';
            if (cleanObj.partitionObj.msPath.empty()) {
                os << "    *** MS does not exist ***\n";
            } else {
                os << "  band            " << getInfo().spectralWindow() << '\n';
                os << "  nchan:          " << getInfo().nchan() << "\n";
                os << "  ncorrelations:  " << getInfo().ncorr() << '\n';
                unsigned int nrbl = getInfo().nbaselines();
                os << "  nbaselines:     " << nrbl << '\n';
                os << "  ntimes:         " << getInfo().ntime() << '\n';  // itsSelMS can contain timeslots that are ignored in process
                os << "  time interval:  " << getInfo().timeInterval() << '\n';
            }
        } 

        void Reorder::finish() { 
            cleanObj.partitionObj.postprocess();
            // std::cout << "Finish" << std::endl;
            getNextStep()->finish(); 
        }


        void Reorder::configureSettings(const common::ParameterSet& parset)
        {
            // cleanObj._settings.divideChannelFrequencies.assign(
            //     parset.getDoubleVector("divideChannelFrequencies", true).begin(),
            //     parset.getDoubleVector("divideChannelFrequencies", true).end());
            // cleanObj._settings.spectralCorrection = parset.getFloatVector("spectralCorrection", true);

            // cleanObj._settings.continuedRun = parset.getBool("continuedRun");
            // cleanObj._settings.parallelReordering = parset.getInt32("parallelReordering");

            // cleanObj._settings.simulatedBaselineNoiseFilename = parset.getString("simulatedBaselineNoiseFilename");
            // cleanObj._settings.atermConfigFilename = parset.getString("atermConfigFilename");
            // cleanObj._settings.ddPsfGridHeight = parset.getInt32("ddPsfGridHeight");
            // cleanObj._settings.ddPsfGridWidth = parset.getInt32("ddPsfGridWidth");
            // cleanObj._settings.subtractModel  = parset.getBool("subtractModel");
            // cleanObj._settings.modelUpdateRequired = parset.getBool("modelUpdateRequired");
            // cleanObj._settings.diagonalSolutions = parset.getBool("diagonalSolutions");
            // cleanObj._settings.gridWithBeam = parset.getBool("gridWithBeam");
            // cleanObj._settings.simulateNoise = parset.getBool("simulateNoise");
            // cleanObj._settings.baselineDependentAveragingInWavelengths = parset.getFloat("baselineDependentAveragingInWavelengths");
            // cleanObj._settings.simulatedNoiseStdDev = parset.getFloat("simulatedNoiseStdDev");
            // cleanObj._settings.deconvolutionMGain = parset.getFloat("deconvolutionMGain");
            
            // cleanObj._settings.hasShift = parset.getBool("hasShift");
            // cleanObj._settings.joinedFrequencyDeconvolution = parset.getBool("joinedFrequencyDeconvolution");
            // cleanObj._settings.joinedPolarizationDeconvolution = parset.getBool("joinedPolarizationDeconvolution");
            // cleanObj._settings.divideChannelsByGaps = parset.getBool("divideChannelsByGaps");
            
            // cleanObj._settings.facetRegionFilename = parset.getString("facetRegionFilename");
            // cleanObj._settings.dataColumnName = parset.getString("dataColumnName");
            // cleanObj._settings.trimmedImageWidth = parset.getInt32("trimmedImageWidth");
            // cleanObj._settings.trimmedImageHeight = parset.getInt32("trimmedImageHeight");
            // cleanObj._settings.deconvolutionIterationCount = parset.getInt32("deconvolutionIterationCount");
            // cleanObj._settings.intervalsOut = parset.getInt32("intervalsOut");
            // cleanObj._settings.endChannel = parset.getInt32("endChannel");
            // cleanObj._settings.startChannel = parset.getInt32("startChannel");
            // cleanObj._settings.channelsOut = parset.getInt32("channelsOut");
            // cleanObj._settings.pixelScaleX = parset.getFloat("pixelScaleX");
            // cleanObj._settings.pixelScaleY = parset.getFloat("pixelScaleY");
            // cleanObj._settings.imagePadding = parset.getFloat("imagePadding");
            // cleanObj._settings.shiftRA = parset.getFloat("shiftRA");
            // cleanObj._settings.shiftDec = parset.getFloat("shiftDec");
            // cleanObj._settings.spectralCorrectionFrequency = parset.getFloat("spectralCorrectionFrequency");
        }

    }  // namespace steps
}