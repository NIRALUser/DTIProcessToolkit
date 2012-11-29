macro(SetTargetDirectories TARGET_NAME)
 if( DEFINED SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY
 AND DEFINED SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY
 AND DEFINED SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY
   )
   set_target_properties(${TARGET_NAME} PROPERTIES
     RUNTIME_OUTPUT_DIRECTORY ${SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY}
     LIBRARY_OUTPUT_DIRECTORY ${SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY}
     ARCHIVE_OUTPUT_DIRECTORY ${SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY}
     )
 endif()
endmacro()

ADD_EXECUTABLE(nrrdvtk nrrdvtk.cxx)
TARGET_LINK_LIBRARIES(nrrdvtk ${ITK_LIBRARIES})

ADD_EXECUTABLE(dwiprocess dwiprocess.cxx)
TARGET_LINK_LIBRARIES(dwiprocess ${ITK_LIBRARIES})

ADD_EXECUTABLE(labeltable labeltable.cxx)
TARGET_LINK_LIBRARIES(labeltable ${ITK_LIBRARIES})

ADD_EXECUTABLE(hfieldtodeformation  hfieldtodeformation.cxx)
TARGET_LINK_LIBRARIES(hfieldtodeformation TensorOperations DTIIO ${ITK_LIBRARIES})

ADD_EXECUTABLE(tio tio.cxx)
TARGET_LINK_LIBRARIES(tio ${ITK_LIBRARIES})

ADD_EXECUTABLE(morph morph.cxx)
TARGET_LINK_LIBRARIES(morph ${ITK_LIBRARIES})

if(BUILD_dwiAtlas)
  ADD_SUBDIRECTORY(dwiAtlas)
  ADD_DEPENDENCIES(dwiAtlas dwiprocess DTIIO)	
endif(BUILD_dwiAtlas)

set(dtiprocess_SRC dtiprocess.cxx)
GENERATECLP(dtiprocess_SRC dtiprocess.xml)
ADD_EXECUTABLE(dtiprocess ${dtiprocess_SRC})
TARGET_LINK_LIBRARIES(dtiprocess TensorOperations DTIIO ${ITK_LIBRARIES})
SetTargetDirectories(dtiprocess)
IF(BUILD_TESTING)
  ADD_LIBRARY( dtiprocessLib STATIC ${dtiprocess_SRC})
  set_target_properties (dtiprocessLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
  TARGET_LINK_LIBRARIES(dtiprocessLib TensorOperations DTIIO ${ITK_LIBRARIES})
ENDIF(BUILD_TESTING)

set(dtiestim_SRC dtiestim.cxx)
GENERATECLP(dtiestim_SRC dtiestim.xml)
ADD_EXECUTABLE(dtiestim ${dtiestim_SRC})
TARGET_LINK_LIBRARIES(dtiestim DWIPrivateLibrary ${ITK_LIBRARIES} cephes)
SetTargetDirectories(dtiestim)
IF(BUILD_TESTING)
  ADD_LIBRARY( dtiestimLib STATIC ${dtiestim_SRC})
  set_target_properties (dtiestimLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
  TARGET_LINK_LIBRARIES(dtiestimLib DWIPrivateLibrary ${ITK_LIBRARIES} cephes)
ENDIF(BUILD_TESTING)


set(scalartransform_SRC scalartransform.cxx)
GENERATECLP(scalartransform_SRC scalartransform.xml)
ADD_EXECUTABLE(scalartransform ${scalartransform_SRC})
TARGET_LINK_LIBRARIES(scalartransform DTIIO ${ITK_LIBRARIES} )
SetTargetDirectories(scalartransform)

set(fiberprocess_SRC fiberprocess.cxx)
GENERATECLP(fiberprocess_SRC fiberprocess.xml)
ADD_EXECUTABLE(fiberprocess ${fiberprocess_SRC})
TARGET_LINK_LIBRARIES(fiberprocess DTIIO ${ITK_LIBRARIES} )
SetTargetDirectories(fiberprocess)

set(dtiaverage_SRC dtiaverage.cxx)
GENERATECLP(dtiaverage_SRC dtiaverage.xml)
ADD_EXECUTABLE(dtiaverage ${dtiaverage_SRC})
TARGET_LINK_LIBRARIES(dtiaverage TensorOperations DTIIO )
SetTargetDirectories(dtiaverage)
IF(BUILD_TESTING)
  ADD_LIBRARY( dtiaverageLib STATIC ${dtiaverage_SRC})
  set_target_properties (dtiaverageLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
  TARGET_LINK_LIBRARIES(dtiaverageLib TensorOperations DTIIO )
ENDIF(BUILD_TESTING)

set(fibertrack_SRC fibertrack.cxx)
GENERATECLP(fibertrack fibertrack.xml)
ADD_EXECUTABLE(fibertrack ${fibertrack_SRC})
TARGET_LINK_LIBRARIES(fibertrack DTIIO TensorOperations ${ITK_LIBRARIES} )
SetTargetDirectories(fibertrack)

set(fiberstats_SRC fiberstats.cxx)
GENERATECLP(fiberstats_SRC fiberstats.xml)
ADD_EXECUTABLE(fiberstats fiberstats.cxx)
TARGET_LINK_LIBRARIES(fiberstats DTIIO )
SetTargetDirectories(fiberstats)

set(maxcurvature_SRC maxcurvature.cxx)
GENERATECLP(maxcurvature_SRC maxcurvature.xml)
ADD_EXECUTABLE(maxcurvature ${maxcurvature_SRC})
TARGET_LINK_LIBRARIES(maxcurvature ${ITK_LIBRARIES} )
SetTargetDirectories(maxcurvature)

INSTALL(TARGETS dtiprocess 
                dtiestim
                scalartransform
                fiberprocess
                dtiaverage
                fibertrack
                fiberstats
                maxcurvature
                 RUNTIME DESTINATION ${INSTALL_RUNTIME_DESTINATION}
                 LIBRARY DESTINATION ${INSTALL_LIBRARY_DESTINATION}
                 ARCHIVE DESTINATION ${INSTALL_ARCHIVE_DESTINATION})

