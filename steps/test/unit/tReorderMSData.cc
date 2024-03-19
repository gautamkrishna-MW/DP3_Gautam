
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <limits>

#include "dp3/base/DPBuffer.h"
#include "dp3/base/DPInfo.h"
#include "dp3/steps/Step.h"
#include "tStepCommon.h"
#include "../../InputStep.h"
#include "../../NullStep.h"
#include "../../ReorderMSData.h"
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::Fields;
using dp3::common::ParameterSet;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(reorder)

// Function to compare two binary files. If the files for comparison are Meta
// file, the flag isMetaFile is set to true.
std::pair<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
compareBinaryFiles(const std::string& referenceFile,
                   const std::string& inputFile, bool isMetaFile) {
  std::pair<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
      filePtrPair;

  std::ifstream refFilePtr(
      referenceFile,
      std::ifstream::ate | std::ifstream::binary);  // open file at the end
  std::ifstream inpFilePtr(
      inputFile,
      std::ifstream::ate | std::ifstream::binary);  // open file at the end

  if (refFilePtr.fail() || inpFilePtr.fail()) {
    std::cout << "Failed to open file" << std::endl;
    return filePtrPair;
  }

  // Resetting file pointers to compare from the beginning.
  refFilePtr.seekg(0);  // rewind
  inpFilePtr.seekg(0);  // rewind

  // Meta file contain path information in the first line of the binary.
  // Since the path representation differs (i.e absolute path vs relative path),
  // we skip the path check and compare only the rest of the meta information.
  if (isMetaFile) {
    refFilePtr.ignore(std::numeric_limits<std::streamsize>::max(),
                      '\n');  // Ignore first line
    inpFilePtr.ignore(std::numeric_limits<std::streamsize>::max(),
                      '\n');  // Ignore first line

    if (refFilePtr.eof() || inpFilePtr.eof()) {
      std::cout << "End of file reached" << std::endl;
      return filePtrPair;
    }
  } else {
    // If file sizes are unequal return false.
    if (refFilePtr.tellg() != inpFilePtr.tellg()) {
      std::cout << "Unequal binary file sizes" << std::endl;
      return filePtrPair;  // different file size
    }
    // Resetting file pointers to compare from the beginning.
    refFilePtr.seekg(0);  // rewind
    inpFilePtr.seekg(0);  // rewind
  }

  std::istreambuf_iterator<char> refFileIterator(refFilePtr);
  std::istreambuf_iterator<char> inpFileIterator(inpFilePtr);

  filePtrPair.first = refFileIterator;
  filePtrPair.second = inpFileIterator;

  return filePtrPair;
}

// Test simple averaging without flagged points.
void test1(std::string msPath, std::string msOutPath) {
  // Create the steps.
  ParameterSet parset;
  parset.add("msin", msPath);
  parset.add("msin.ntimes", "2");

  if (msOutPath.empty()) {
    std::string tmpPath(msPath);
    tmpPath.resize(tmpPath.size() - 3);
    parset.add("msout", tmpPath + "_output.ms");
  } else
    parset.add("msout", msOutPath);

  // Creating steps
  std::shared_ptr<dp3::steps::InputStep> input_step =
      std::move(dp3::steps::InputStep::CreateReader(parset));
  auto reorderStep = std::make_shared<dp3::steps::Reorder>(parset, "");
  auto nullStep = std::make_shared<dp3::steps::NullStep>();

  // Assiging the fields to read
  input_step->setFieldsToRead(Step::kDataField | Step::kFlagsField |
                              Step::kWeightsField | Step::kUvwField);
  dp3::steps::test::Execute({input_step, reorderStep, nullStep});
}

BOOST_AUTO_TEST_CASE(testReorder1) {
  test1("../resources/midbands/reorderTestMS.ms", "./testOut.ms");

  // Default stream buffer iterator has a end-of-stream state, so using this to
  // detect end of the stream.
  std::istreambuf_iterator<char> endOfStream;

  // Data File check
  std::pair<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
      filePtrPair_vis;
  filePtrPair_vis =
      compareBinaryFiles("../resources/midbands/refTmps/reorderTestVis.tmp",
                         "./testOut.ms-part0000-I-b0.tmp", false);
  BOOST_REQUIRE_EQUAL_COLLECTIONS(filePtrPair_vis.first, endOfStream,
                                  filePtrPair_vis.second, endOfStream);

  // Weight File check
  std::pair<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
      filePtrPair_weight;
  filePtrPair_weight =
      compareBinaryFiles("../resources/midbands/refTmps/reorderTestWeights.tmp",
                         "./testOut.ms-part0000-I-b0-w.tmp", false);
  BOOST_REQUIRE_EQUAL_COLLECTIONS(filePtrPair_weight.first, endOfStream,
                                  filePtrPair_weight.second, endOfStream);

  // Meta File check
  std::pair<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
      filePtrPair_meta;
  filePtrPair_meta =
      compareBinaryFiles("../resources/midbands/refTmps/reorderTestMeta.tmp",
                         "./testOut.ms-spw0-parted-meta.tmp", true);
  BOOST_REQUIRE_EQUAL_COLLECTIONS(filePtrPair_meta.first, endOfStream,
                                  filePtrPair_meta.second, endOfStream);
}

BOOST_AUTO_TEST_SUITE_END()