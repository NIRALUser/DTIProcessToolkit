add_executable(nrrdvtk nrrdvtk.cxx)
target_link_libraries(nrrdvtk ${ITK_LIBRARIES})

add_executable(dwiprocess dwiprocess.cxx)
target_link_libraries(dwiprocess ${ITK_LIBRARIES})

add_executable(labeltable labeltable.cxx)
target_link_libraries(labeltable ${ITK_LIBRARIES})

add_executable(hfieldtodeformation  hfieldtodeformation.cxx)
target_link_libraries(hfieldtodeformation TensorOperations DTIIO ${ITK_LIBRARIES})

add_executable(tio tio.cxx)
target_link_libraries(tio ${ITK_LIBRARIES})

add_executable(morph morph.cxx)
target_link_libraries(morph ${ITK_LIBRARIES})

if(BUILD_dwiAtlas)
  add_subdirectory(dwiAtlas)
  add_dependencies(dwiAtlas dwiprocess DTIIO)
endif(BUILD_dwiAtlas)


set(dtiprocess_SRC dtiprocess.cxx)
GENERATECLP(dtiprocess_SRC dtiprocess.xml)
add_executable(dtiprocess ${dtiprocess_SRC})
target_link_libraries(dtiprocess TensorOperations DTIIO ${ITK_LIBRARIES})

set(dtiestim_SRC dtiestim.cxx)
GENERATECLP(dtiestim_SRC dtiestim.xml)
add_executable(dtiestim ${dtiestim_SRC})
target_link_libraries(dtiestim ${ITK_LIBRARIES} cephes)

set(scalartransform_SRC scalartransform.cxx)
GENERATECLP(scalartransform_SRC scalartransform.xml)
add_executable(scalartransform ${scalartransform_SRC})
target_link_libraries(scalartransform DTIIO ${ITK_LIBRARIES} )

set(fiberprocess_SRC fiberprocess.cxx)
GENERATECLP(fiberprocess_SRC fiberprocess.xml)
add_executable(fiberprocess ${fiberprocess_SRC})
target_link_libraries(fiberprocess DTIIO ${ITK_LIBRARIES} ${ITK_LIBRARIES} )

set(dtiaverage_SRC dtiaverage.cxx)
GENERATECLP(dtiaverage_SRC dtiaverage.xml)
add_executable(dtiaverage ${dtiaverage_SRC})
target_link_libraries(dtiaverage TensorOperations DTIIO )

set(fibertrack_SRC fibertrack.cxx)
GENERATECLP(fibertrack fibertrack.xml)
add_executable(fibertrack ${fibertrack_SRC})
target_link_libraries(fibertrack DTIIO TensorOperations ${ITK_LIBRARIES} ${ITK_LIBRARIES})

set(fiberstats_SRC fiberstats.cxx)
GENERATECLP(fiberstats_SRC fiberstats.xml)
add_executable(fiberstats fiberstats.cxx)
target_link_libraries(fiberstats DTIIO )

set(maxcurvature_SRC maxcurvature.cxx)
GENERATECLP(maxcurvature_SRC maxcurvature.xml)
add_executable(maxcurvature ${maxcurvature_SRC})
target_link_libraries(maxcurvature ${ITK_LIBRARIES} ${ITK_LIBRARIES} )


install(TARGETS dtiprocess
                dtiestim
                scalartransform
                #                fiberprocess
                dtiaverage
                fibertrack
                fiberstats
                maxcurvature
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib/static)

