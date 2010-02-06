
# Unit tests
SET(TestSuite_SRCS
BesselTest.cxx
ScalarDeformationApplyTest.cxx
WarpFiberTest.cxx
WarpTensorTest.cxx
FiberIOTest.cxx
TensorEstimateTest.cxx
PrincipalEigenvectorFunctionTest.cxx
ReadITKAffineTest.cxx
)

ADD_EXECUTABLE(TestSuite TestSuite.cxx ${TestSuite_SRCS})
TARGET_LINK_LIBRARIES(TestSuite ITKNumerics ITKCommon ITKIO DWIPrivateLibrary)


SET(TESTING_DATA_DIR ${CMAKE_SOURCE_DIR}/Testing/Data)
SET(TESTING_RESULTS_DIR ${CMAKE_BINARY_DIR}/Testing/Data)

FILE(MAKE_DIRECTORY ${TESTING_RESULTS_DIR})

# Unit tests
SET(TEST_SUITE ${CMAKE_BINARY_DIR}/bin/TestSuite)

ADD_TEST(BesselTest ${TEST_SUITE} BesselTest)
ADD_TEST(LLSTensorEstimateTest ${TEST_SUITE} LLSTensorEstimateTest)
ADD_TEST(WLSTensorEstimateTest ${TEST_SUITE} WLSTensorEstimateTest)
ADD_TEST(NLSTensorEstimateTest ${TEST_SUITE} NLSTensorEstimateTest)
ADD_TEST(PrincipalEigenvectorFunctionTest ${TEST_SUITE} PrincipalEigenvectorFunctionTest)
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

ADD_TEST(ReadITKAffineTest ${TEST_SUITE} ReadITKAffineTest
                           ${TESTING_DATA_DIR}/itkAffineTransformFloat.txt
                           ${TESTING_DATA_DIR}/itkAffineTransformDouble.txt)

# Application Command Line tests
ADD_TEST(TransformDWITest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}dwiprocess
                          ${TESTING_DATA_DIR}/testdwi.nrrd
                          ${TESTING_RESULTS_DIR}/testdwi-rot.nrrd
                          ${TESTING_DATA_DIR}/test.dof)

ADD_TEST(TensorEstimateVolumeTest
         ${TEST_SUITE}
         --compare
         ${TESTING_DATA_DIR}/tensors.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors.nrrd
         CommandLineTest
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
         ${TESTING_DATA_DIR}/testdwi.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors.nrrd -t 200)

ADD_TEST(TensorEstimateNEXTest
         ${TEST_SUITE}
         --compare
         ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd
         CommandLineTest
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
         ${TESTING_DATA_DIR}/testdwirepbaseline.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd -t 200 &&
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
         ${TESTING_DATA_DIR}/testdwirepbaselinenex.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd -t 200)

ADD_TEST(TensorEstimateVolumeTest
         ${TEST_SUITE}
         --compare
         ${TESTING_DATA_DIR}/tensors.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors.nrrd
         CommandLineTest
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
         ${TESTING_DATA_DIR}/testdwi.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors.nrrd -t 200)

ADD_TEST(TensorEstimateNEXTest
         ${TEST_SUITE}
         --compare
         ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd
         CommandLineTest
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
         ${TESTING_DATA_DIR}/testdwirepbaseline.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd -t 200 &&
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
         ${TESTING_DATA_DIR}/testdwirepbaselinenex.nrrd
         ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd -t 200)

ADD_TEST(CompareFATest
         ${TEST_SUITE}
         --compare ${TESTING_DATA_DIR}/fa.nrrd
                   ${TESTING_RESULTS_DIR}/result-fa.nrrd
         CommandLineTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd -f ${TESTING_RESULTS_DIR}/result-fa.nrrd)

# ADD_TEST(CompareIntFAPythonTest
#          ${TEST_SUITE}
#          --compare ${TESTING_DATA_DIR}/fa.nrrd
#                    ${TESTING_RESULTS_DIR}/presult-fa.nrrd
#          CommandLineTest PYTHONPATH=${CMAKE_BINARY_DIR}/Wrapping python ${CMAKE_SOURCE_DIR}/Testing/test.py ${TESTING_DATA_DIR}/tensors.nrrd  ${TESTING_RESULTS_DIR}/presult-fa.nrrd int)

# ADD_TEST(CompareFloatFAPythonTest
#          ${TEST_SUITE}
#          --compare ${TESTING_DATA_DIR}/floatfa.nrrd
#                    ${TESTING_RESULTS_DIR}/presult-float-fa.nrrd
#          CommandLineTest PYTHONPATH=${CMAKE_BINARY_DIR}/Wrapping python ${CMAKE_SOURCE_DIR}/Testing/test.py ${TESTING_DATA_DIR}/tensors.nrrd  ${TESTING_RESULTS_DIR}/presult-float-fa.nrrd float)

ADD_TEST(CompareCFATest
         ${TEST_SUITE}
         --compare ${TESTING_DATA_DIR}/cfa.nrrd
                   ${TESTING_RESULTS_DIR}/result-cfa.nrrd
         CommandLineTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd -c ${TESTING_RESULTS_DIR}/result-cfa.nrrd)

ADD_TEST(CompareMDTest
         ${TEST_SUITE}
         --compare ${TESTING_DATA_DIR}/md.nrrd
                   ${TESTING_RESULTS_DIR}/result-md.nrrd
         CommandLineTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd -m ${TESTING_RESULTS_DIR}/result-md.nrrd)

ADD_TEST(CompareLbd1Test
         ${TEST_SUITE}
         --compare ${TESTING_DATA_DIR}/l1.nrrd
                   ${TESTING_RESULTS_DIR}/result-l1.nrrd
         CommandLineTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd --lambda1-output ${TESTING_RESULTS_DIR}/result-l1.nrrd)

ADD_TEST(CompareLbd2Test
         ${TEST_SUITE}
         --compare ${TESTING_DATA_DIR}/l2.nrrd
                   ${TESTING_RESULTS_DIR}/result-l2.nrrd
         CommandLineTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd --lambda2-output ${TESTING_RESULTS_DIR}/result-l2.nrrd)

ADD_TEST(CompareLbd3Test
         ${TEST_SUITE}
         --compare ${TESTING_DATA_DIR}/l3.nrrd
                   ${TESTING_RESULTS_DIR}/result-l3.nrrd
         CommandLineTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess ${TESTING_DATA_DIR}/tensors.nrrd --lambda3-output ${TESTING_RESULTS_DIR}/result-l3.nrrd)

