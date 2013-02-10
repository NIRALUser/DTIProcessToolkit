
# Unit tests
SET(DTIProcessTestHarness_SRCS
BesselTest.cxx
ScalarDeformationApplyTest.cxx
WarpFiberTest.cxx
WarpTensorTest.cxx
FiberIOTest.cxx
TensorEstimateTest.cxx
PrincipalEigenvectorFunctionTest.cxx
ReadITKAffineTest.cxx
)

ADD_EXECUTABLE(DTIProcessTestHarness DTIProcessTestHarness.cxx ${DTIProcessTestHarness_SRCS})
TARGET_LINK_LIBRARIES(DTIProcessTestHarness ${ITK_LIBRARIES} cephes)

#SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})

SET(TESTING_DATA_DIR ${CMAKE_SOURCE_DIR}/Testing/Data)
SET(TESTING_RESULTS_DIR ${CMAKE_BINARY_DIR}/Testing/Data)

FILE(MAKE_DIRECTORY ${TESTING_RESULTS_DIR})

# Unit tests

add_test(NAME BesselTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> BesselTest)
add_test(NAME LLSTensorEstimateTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> LLSTensorEstimateTest)
add_test(NAME WLSTensorEstimateTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> WLSTensorEstimateTest)
add_test(NAME NLSTensorEstimateTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> NLSTensorEstimateTest)
add_test(NAME PrincipalEigenvectorFunctionTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> PrincipalEigenvectorFunctionTest)
#add_test(NAME WarpFiberTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> WarpFiberTest)
#add_test(NAME WarpTensortest COMMAND $<TARGET_FILE:DTIProcessTestHarness> WarpTensorTest)

add_test(NAME ScalarDeformationApplyTest COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_RESULTS_DIR}/transformfa.nrrd
  ${TESTING_RESULTS_DIR}/warpfa.nrrd
  ScalarDeformationApplyTest
  ${TESTING_DATA_DIR}/fa.nrrd
  ${TESTING_DATA_DIR}/test.dof
  ${TESTING_RESULTS_DIR}/transformfa.nrrd
  ${TESTING_RESULTS_DIR}/warpfa.nrrd
  )

add_test(NAME ReadITKAffineTest COMMAND $<TARGET_FILE:DTIProcessTestHarness> ReadITKAffineTest
  ${TESTING_DATA_DIR}/itkAffineTransformFloat.txt
  ${TESTING_DATA_DIR}/itkAffineTransformDouble.txt)

# Application Command Line tests
add_test(NAME TransformDWITest COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dwiprocess
  ${TESTING_DATA_DIR}/testdwi.nrrd
  ${TESTING_RESULTS_DIR}/testdwi-rot.nrrd
  ${TESTING_DATA_DIR}/test.dof)

add_test(NAME TensorEstimateVolumeTest1 COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/tensors.nrrd
  ${TESTING_RESULTS_DIR}/result-tensors.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
  --dwi_image ${TESTING_DATA_DIR}/testdwi.nrrd
  --tensor_output ${TESTING_RESULTS_DIR}/result-tensors.nrrd
  --threshold 200)

add_test(NAME TensorEstimateNEXTest1 COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd
  ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
  --dwi_image ${TESTING_DATA_DIR}/testdwirepbaseline.nrrd
  --tensor_output ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd
  --threshold 200
  &&
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
  --dwi_image ${TESTING_DATA_DIR}/testdwirepbaselinenex.nrrd
  --tensor_output ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd
  --threshold 200)

add_test(NAME TensorEstimateVolumeTest2 COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/tensors.nrrd
  ${TESTING_RESULTS_DIR}/result-tensors.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
  --dwi_image ${TESTING_DATA_DIR}/testdwi.nrrd
  --tensor_output ${TESTING_RESULTS_DIR}/result-tensors.nrrd
  -t 200)

add_test(NAME TensorEstimateNEXTest2 COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd
  ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
  --dwi_image ${TESTING_DATA_DIR}/testdwirepbaseline.nrrd
  --tensor_output ${TESTING_RESULTS_DIR}/result-tensors-rep.nrrd
  -t 200
  &&
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiestim
  --dwi_image ${TESTING_DATA_DIR}/testdwirepbaselinenex.nrrd
  --tensor_output ${TESTING_RESULTS_DIR}/result-tensors-nex.nrrd
  -t 200)

add_test(NAME CompareFATest COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/fa.nrrd
  ${TESTING_RESULTS_DIR}/result-fa.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess
  --dti_image ${TESTING_DATA_DIR}/tensors.nrrd
  -f ${TESTING_RESULTS_DIR}/result-fa.nrrd)

# add_test(NAME CompareIntFAPythonTest COMMAND $<TARGET_FILE:DTIProcessTestHarness>
#   --compare ${TESTING_DATA_DIR}/fa.nrrd
#   ${TESTING_RESULTS_DIR}/presult-fa.nrrd
#   CommandLineTest PYTHONPATH=${CMAKE_BINARY_DIR}/Wrapping python ${CMAKE_SOURCE_DIR}/Testing/test.py ${TESTING_DATA_DIR}/tensors.nrrd  ${TESTING_RESULTS_DIR}/presult-fa.nrrd int)

# add_test(NAME CompareFloatFAPythonTest COMMAND $<TARGET_FILE:DTIProcessTestHarness>
#   --compare ${TESTING_DATA_DIR}/floatfa.nrrd
#   ${TESTING_RESULTS_DIR}/presult-float-fa.nrrd
#   CommandLineTest PYTHONPATH=${CMAKE_BINARY_DIR}/Wrapping python ${CMAKE_SOURCE_DIR}/Testing/test.py ${TESTING_DATA_DIR}/tensors.nrrd  ${TESTING_RESULTS_DIR}/presult-float-fa.nrrd float)

add_test(NAME CompareCFATest COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/cfa.nrrd
  ${TESTING_RESULTS_DIR}/result-cfa.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess
  --dti_image ${TESTING_DATA_DIR}/tensors.nrrd
  --color_fa_output ${TESTING_RESULTS_DIR}/result-cfa.nrrd)

add_test(NAME CompareMDTest COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/md.nrrd
  ${TESTING_RESULTS_DIR}/result-md.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess --dti_image
  ${TESTING_DATA_DIR}/tensors.nrrd -m
  ${TESTING_RESULTS_DIR}/result-md.nrrd)

add_test(NAME CompareLbd1Test COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/l1.nrrd
  ${TESTING_RESULTS_DIR}/result-l1.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess --dti_image
  ${TESTING_DATA_DIR}/tensors.nrrd --lambda1_output
  ${TESTING_RESULTS_DIR}/result-l1.nrrd)

add_test(NAME CompareLbd2Test COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/l2.nrrd
  ${TESTING_RESULTS_DIR}/result-l2.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess --dti_image
  ${TESTING_DATA_DIR}/tensors.nrrd --lambda2_output
  ${TESTING_RESULTS_DIR}/result-l2.nrrd)

add_test(NAME CompareLbd3Test COMMAND $<TARGET_FILE:DTIProcessTestHarness>
  --compare
  ${TESTING_DATA_DIR}/l3.nrrd
  ${TESTING_RESULTS_DIR}/result-l3.nrrd
  CommandLineTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dtiprocess --dti_image
  ${TESTING_DATA_DIR}/tensors.nrrd --lambda3_output
  ${TESTING_RESULTS_DIR}/result-l3.nrrd)

