
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
void compareBinaryFiles(const std::string& reference_filename,
                        const std::string& input_filename, bool isMetaFile) {
  std::pair<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>>
      iterator_pair;

  std::ifstream reference_file_ptr(reference_filename);
  std::ifstream input_file_ptr(input_filename);

  if (reference_file_ptr.fail() || input_file_ptr.fail()) {
    std::cout << "Failed to open file" << std::endl;
    BOOST_CHECK(false);
  }

  // Meta files contain path information in the first line of the binary.
  // Since the path representation differs (i.e absolute path vs relative path),
  // we skip the path check and compare only the rest of the meta information.
  if (isMetaFile) {
    reference_file_ptr.ignore(std::numeric_limits<std::streamsize>::max(),
                              '\n');
    input_file_ptr.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if (reference_file_ptr.eof() || input_file_ptr.eof()) {
      std::cout << "End of file reached" << std::endl;
      BOOST_CHECK(false);
    }
  } else {
    // If file sizes are unequal return false.
    if (reference_file_ptr.tellg() != input_file_ptr.tellg()) {
      std::cout << "Unequal binary file sizes" << std::endl;
      BOOST_CHECK(false);
    }

    // Resetting file pointers to compare from the beginning.
    reference_file_ptr.seekg(0);
    input_file_ptr.seekg(0);
  }

  // Setting the file stream to not to skip whitespaces.
  reference_file_ptr >> std::noskipws;
  input_file_ptr >> std::noskipws;
  std::istreambuf_iterator<char> reference_file_iterator(reference_file_ptr);
  std::istreambuf_iterator<char> input_file_iterator(input_file_ptr);

  std::istreambuf_iterator<char> end_of_stream;
  BOOST_REQUIRE_EQUAL_COLLECTIONS(reference_file_iterator, end_of_stream,
                                  input_file_iterator, end_of_stream);
}

// Test simple averaging without flagged points.
void ReorderTest(std::string ms_input_path, std::string ms_output_path) {
  // Create the steps.
  ParameterSet parset;
  parset.add("msin", ms_input_path);
  parset.add("msin.ntimes", "2");

  if (ms_output_path.empty()) {
    std::string temp_path(ms_input_path);
    temp_path.resize(temp_path.size() - 3);
    parset.add("msout", temp_path + "_output.ms");
  } else
    parset.add("msout", ms_output_path);

  // Creating steps
  std::shared_ptr<dp3::steps::InputStep> input_step =
      std::move(dp3::steps::InputStep::CreateReader(parset));
  auto reorder_step = std::make_shared<dp3::steps::Reorder>(parset, "");
  auto null_step = std::make_shared<dp3::steps::NullStep>();

  // Assiging the fields to read
  input_step->setFieldsToRead(Step::kDataField | Step::kFlagsField |
                              Step::kWeightsField | Step::kUvwField);
  dp3::steps::test::Execute({input_step, reorder_step, null_step});
}

BOOST_AUTO_TEST_CASE(testReorder1) {
  ReorderTest("../resources/midbands/reorderTestMS.ms", "./testOut.ms");

  // Data File check
  compareBinaryFiles("../resources/midbands/refTmps/reorderTestVis.tmp",
                     "./testOut.ms-part0000-I-b0.tmp", false);

  // Weight File check
  compareBinaryFiles("../resources/midbands/refTmps/reorderTestWeights.tmp",
                     "./testOut.ms-part0000-I-b0-w.tmp", false);

  // Meta File check
  compareBinaryFiles("../resources/midbands/refTmps/reorderTestMeta.tmp",
                     "./testOut.ms-spw0-parted-meta.tmp", true);
}

BOOST_AUTO_TEST_SUITE_END()