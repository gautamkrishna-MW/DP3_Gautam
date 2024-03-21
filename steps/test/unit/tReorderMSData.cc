
#include <algorithm>
#include <fstream>
#include <limits>

#include <boost/test/unit_test.hpp>

#include "tStepCommon.h"
#include "dp3/steps/Step.h"
#include "../../InputStep.h"
#include "../../NullStep.h"
#include "../../ReorderMSData.h"
#include "../../../common/ParameterSet.h"

using dp3::common::Fields;
using dp3::common::ParameterSet;
using dp3::steps::Step;

BOOST_AUTO_TEST_SUITE(reorder)

// Function to compare two binary files. If the files for comparison are Meta
// file, the flag is_meta_file is set to true.
void compareBinaryFiles(const std::string& reference_filename,
                        const std::string& input_filename, bool is_meta_file) {
  
  std::ifstream reference_file(reference_filename);
  std::ifstream input_file(input_filename);

  if (reference_file.fail() || input_file.fail()) {
    BOOST_CHECK(false);
  }

  // Meta files contain path information in the first line of the binary.
  // Since the path representation differs (i.e absolute path vs relative path),
  // we skip the path check and compare only the rest of the meta information.
  if (is_meta_file) {
    reference_file.ignore(std::numeric_limits<std::streamsize>::max(),
                              '\n');
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    if (reference_file.eof() || input_file.eof()) {
      BOOST_CHECK(false);
    }
  } else {
    // If file sizes are unequal return false.
    if (reference_file.tellg() != input_file.tellg()) {
      BOOST_CHECK(false);
    }

    // Resetting file pointers to compare from the beginning.
    reference_file.seekg(0);
    input_file.seekg(0);
  }

  // Creating file stream iterators and check if they are equal
  std::istreambuf_iterator<char> reference_file_iterator(reference_file);
  std::istreambuf_iterator<char> input_file_iterator(input_file);

  std::istreambuf_iterator<char> end_of_stream;
  BOOST_REQUIRE_EQUAL_COLLECTIONS(reference_file_iterator, end_of_stream,
                                  input_file_iterator, end_of_stream);
}

// Test reorder functionality with default settings
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
  std::shared_ptr<dp3::steps::Reorder> reorder_step =
      std::make_shared<dp3::steps::Reorder>(parset, "");
  std::shared_ptr<dp3::steps::NullStep> null_step =
      std::make_shared<dp3::steps::NullStep>();

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