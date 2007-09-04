#ADD_TEST(bessel)
SET(TESTING_DATA_DIR ${CMAKE_SOURCE_DIR}/Testing/Data)
SET(TESTING_RESULTS_DIR ${CMAKE_BINARY_DIR}/Testing/Data)

FILE(MAKE_DIRECTORY ${TESTING_RESULTS_DIR})

# Unit tests
SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
SET(TEST_SUITE ${CXX_TEST_PATH}/TestSuite)

ADD_TEST(BesselTest ${TEST_SUITE} BesselTest)
ADD_TEST(TensorEstimateTest ${TEST_SUITE} TensorEstimateTest)
#ADD_TEST(WarpFiberTest ${TEST_SUITE} WarpFiberTest)
#ADD_TEST(WarpTensortest ${TEST_SUITE} WarpTensorTest)

ADD_TEST(ScalarDeformationApplyTest ${TEST_SUITE}
         --compare
           ${TESTING_RESULTS_DIR}/transformfa.nrrd
           ${TESTING_RESULTS_DIR}/warpfa.nrrd
         ScalarDeformationApplyTest
           ${TESTING_DATA_DIR}/fa.nrrd
           ${TESTING_DATA_DIR}/test.dof
           ${TESTING_RESULTS_DIR}/transformfa.nrrd
           ${TESTING_RESULTS_DIR}/warpfa.nrrd
                     )

# Application Command Line tests
ADD_TEST(TransformDWITest ${CMAKE_BINARY_DIR}/bin/dwiprocess ${TESTING_DATA_DIR}/testdwi.nrrd ${TESTING_RESULTS_DIR}/testdwi-rot.nrrd ${TESTING_DATA_DIR}/test.dof)

ADD_TEST(CompareFATest 
         ${TEST_SUITE} 
         --compare ${TESTING_DATA_DIR}/fa.nrrd
                   ${TESTING_RESULTS_DIR}/result-fa.nrrd
         CommandLineTest ${EXECUTABLE_OUTPUT_PATH}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd -f ${TESTING_RESULTS_DIR}/result-fa.nrrd)

ADD_TEST(CompareCFATest 
         ${TEST_SUITE} 
         --compare ${TESTING_DATA_DIR}/cfa.nrrd
                   ${TESTING_RESULTS_DIR}/result-cfa.nrrd
         CommandLineTest ${EXECUTABLE_OUTPUT_PATH}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd -c ${TESTING_RESULTS_DIR}/result-cfa.nrrd)

ADD_TEST(CompareMDTest 
         ${TEST_SUITE} 
         --compare ${TESTING_DATA_DIR}/md.nrrd
                   ${TESTING_RESULTS_DIR}/result-md.nrrd
         CommandLineTest ${EXECUTABLE_OUTPUT_PATH}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd -m ${TESTING_RESULTS_DIR}/result-md.nrrd)


# Unit tests
SET(TestSuite_SRCS
BesselTest.cxx
ScalarDeformationApplyTest.cxx
#WarpFiberTest.cxx
#WarpTensorTest.cxx
FiberIOTest.cxx
TensorEstimateTest.cxx
)

ADD_EXECUTABLE(TestSuite TestSuite.cxx ${TestSuite_SRCS})
TARGET_LINK_LIBRARIES(TestSuite ITKNumerics ITKCommon ITKIO DWIPrivateLibrary)
