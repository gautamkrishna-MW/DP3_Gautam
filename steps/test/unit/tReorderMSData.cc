
#include <fstream>
#include <algorithm>

#include "tStepCommon.h"
#include "../../ReorderMSData.h"
#include <dp3/base/DPBuffer.h>
#include <dp3/base/DPInfo.h>
#include "../../../common/ParameterSet.h"
#include "../../../common/StringTools.h"

using dp3::base::DPBuffer;
using dp3::base::DPInfo;
using dp3::common::ParameterSet;
using dp3::steps::Step;

bool compareBinaries(const std::string& filename1, const std::string& filename2)
{
    std::ifstream file1(filename1, std::ifstream::ate | std::ifstream::binary); //open file at the end
    std::ifstream file2(filename2, std::ifstream::ate | std::ifstream::binary); //open file at the end
    const std::ifstream::pos_type fileSize = file1.tellg();

    if (fileSize != file2.tellg()) {
        return false; //different file size
    }

    file1.seekg(0); //rewind
    file2.seekg(0); //rewind

    std::istreambuf_iterator<char> begin1(file1);
    std::istreambuf_iterator<char> begin2(file2);

    return std::equal(begin1,std::istreambuf_iterator<char>(),begin2); //Second argument is end-of-range iterator
}

// Test simple averaging without flagged points.
void test1(std::string msPath, std::string msOutPath) {

    // Create the steps.
    ParameterSet parset;
    parset.add("msin", msPath);

    if (msOutPath.empty())
    {
        std:string tmpPath(msPath);
        tmpPath.resize(tmpPath.size()-3);
        parset.add("msout", tmpPath + "_output.ms");
    }
    else
        parset.add("msout", msOutPath);

    auto reorderStep = std::make_shared<dp3::steps::Reorder>(parset, "");
    dp3::steps::test::Execute({reorderStep});
}

BOOST_AUTO_TEST_CASE(testReorder1) { 

    test1("../resources/midbands/midbands_averaged.ms", "./testOut.ms");
    
    // Data File check
    BOOST_CHECK(compareBinaries("../resources/midbands/midbands_averaged.ms-part0000-I-b0.tmp", 
        "./midbands_averaged.ms-part0000-I-b0.tmp"));

    // Weight File check
    BOOST_CHECK(compareBinaries("../resources/midbands_averaged.ms-part0000-I-b0-w.tmp", 
        "./midbands_averaged.ms-part0000-I-b0-w.tmp"));

    // Meta File check
    // BOOST_CHECK(compareBinaries("../resources/midbands_averaged.ms-spw0-parted-meta.tmp", 
    //     "./midbands_averaged.ms-spw0-parted-meta.tmp"));

    // Clean-up
    std::remove("./midbands_averaged.ms-part0000-I-b0.tmp");
    std::remove("./midbands_averaged.ms-part0000-I-b0-w.tmp");
    // std::remove("./midbands_averaged.ms-spw0-parted-meta.tmp");
}