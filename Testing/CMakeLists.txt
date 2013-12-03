#
# Test executable
#



# Add test of the homemade "round()" function
add_executable(TestHomemadeRoundFunction TestHomemadeRoundFunction.cxx)
target_link_libraries(TestHomemadeRoundFunction ${ITK_LIBRARIES})
add_test(NAME TestHomemadeRoundFunction COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:TestHomemadeRoundFunction> )

set(SOURCE_DIRECTORY ${DTIProcess_SOURCE_DIR}/Data/ )
set(TEMP_DIR ${DTIProcess_BINARY_DIR}/Testing/Temporary )

######################################
# DTIEstim tests
######################################

set( CLP dtiestim )
set( ${CLP}_tmp_dir ${TEMP_DIR}/${CLP} )
set( ${CLP}_source_dir ${SOURCE_DIRECTORY}/${CLP} )
file(MAKE_DIRECTORY  ${${CLP}_tmp_dir} )

set(input ${${CLP}_source_dir}/Input/dwi.nrrd )

add_executable(${CLP}Test ImageCompareTest.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)
set(trash ${${CLP}_tmp_dir}/dti_trash.nrrd )

#B0
set(output ${${CLP}_tmp_dir}/b0.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/b0.nrrd )
add_test(NAME ${CLP}B0Test COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0
  ModuleEntryPoint
    --B0 ${output}
    --dwi_image ${input}
    --tensor_output ${trash}
  )

#DWI wls
set(output ${${CLP}_tmp_dir}/dti_wls.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/dti_WLS.nrrd )
add_test(NAME ${CLP}DTI_WLS_Test COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0.000000000001
  ModuleEntryPoint
    --tensor_output ${output}
    --dwi_image ${input}
    -m wls
    --threshold 0
  )

set(output ${${CLP}_tmp_dir}/idwi.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/IDWI.nrrd )
add_test(NAME ${CLP}IDWITest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0.0000000001
  ModuleEntryPoint
    --idwi ${output}
    --dwi_image ${input}
    --tensor_output ${trash}
  )


######################################
# DTIAverage tests
######################################
set( CLP dtiaverage )
set( ${CLP}_tmp_dir ${TEMP_DIR}/${CLP} )
set( ${CLP}_source_dir ${SOURCE_DIRECTORY}/${CLP} )
file(MAKE_DIRECTORY  ${${CLP}_tmp_dir} )

set(input1 ${${CLP}_source_dir}/Input/dti.nrrd )
set(input2 ${${CLP}_source_dir}/Input/dti2.nrrd )

add_executable(${CLP}Test ImageCompareTest.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)


set(output ${${CLP}_tmp_dir}/dti_avg.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/dti_average.nrrd )
add_test(NAME ${CLP}Test1 COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0.000000000000000000001
  ModuleEntryPoint
    --tensor_output ${output}
    --inputs ${input1}
    --inputs ${input2}
  )

######################################
# DTIProcess tests
######################################


set( CLP dtiprocess )
set( ${CLP}_tmp_dir ${TEMP_DIR}/${CLP} )
set( ${CLP}_source_dir ${SOURCE_DIRECTORY}/${CLP} )
file(MAKE_DIRECTORY  ${${CLP}_tmp_dir} )

set(input ${${CLP}_source_dir}/Input/dti.nrrd )

add_executable(${CLP}Test ImageCompareTest.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)

#FA
set(output ${${CLP}_tmp_dir}/fa.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/FA.nrrd )
add_test(NAME ${CLP}FATest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0
  ModuleEntryPoint
    -f ${output}
    --dti_image ${input}
  )

#MD
set(output ${${CLP}_tmp_dir}/MD.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/MD.nrrd )
add_test(NAME ${CLP}MDTest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0
  ModuleEntryPoint
    -m ${output}
    --dti_image ${input}
  )

#AD
set(output ${${CLP}_tmp_dir}/AD.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/AD.nrrd )
add_test(NAME ${CLP}ADTest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0
  ModuleEntryPoint
    --lambda1_output ${output}
    --dti_image ${input}
  )

#RD
set(output ${${CLP}_tmp_dir}/RD.nrrd )
set(baseline ${${CLP}_source_dir}/Baseline/RD.nrrd )
add_test(NAME ${CLP}RDTest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}Test>
  --compare
    ${baseline}
    ${output}
  --compareIntensityTolerance 0
  ModuleEntryPoint
    --RD_output ${output}
    --dti_image ${input}
  )

