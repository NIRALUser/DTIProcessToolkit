ADD_EXECUTABLE(nrrdvtk nrrdvtk.cxx)
TARGET_LINK_LIBRARIES(nrrdvtk ${DTIProcess_ITK_LIBRARIES})
export(TARGETS nrrdvtk FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)

ADD_EXECUTABLE(dwiprocess dwiprocess.cxx)
TARGET_LINK_LIBRARIES(dwiprocess ${DTIProcess_ITK_LIBRARIES})
export(TARGETS dwiprocess APPEND FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)

ADD_EXECUTABLE(labeltable labeltable.cxx)
TARGET_LINK_LIBRARIES(labeltable ${DTIProcess_ITK_LIBRARIES})
export(TARGETS labeltable APPEND FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)

ADD_EXECUTABLE(hfieldtodeformation  hfieldtodeformation.cxx)
TARGET_LINK_LIBRARIES(hfieldtodeformation TensorOperations DTIIO ${DTIProcess_ITK_LIBRARIES})
export(TARGETS hfieldtodeformation APPEND FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)

ADD_EXECUTABLE(tio tio.cxx)
TARGET_LINK_LIBRARIES(tio ${DTIProcess_ITK_LIBRARIES})
export(TARGETS tio APPEND FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)

ADD_EXECUTABLE(morph morph.cxx)
TARGET_LINK_LIBRARIES(morph ${DTIProcess_ITK_LIBRARIES})
export(TARGETS morph APPEND FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)

# When built as a Slicer extension, this will ensure that the CLI executables
# can resolve the path of the ITK/VTK/... libraries provided by Slicer.
#
# It is usually done in UseSlicer.cmake included after calling
# find_package(Slicer)/include(${Slicer_USE_FILE}). But since in this project
# these calls are only done in SuperBuild.cmake, we explicily set the executable
# and library link flags here.
#
# Note: These flags are currently used independently of the value of
#       DTIProcess_BUILD_SLICER_EXTENSION. If this is a problem, the approach
#       should be revisited so that the parameter DTIProcess_BUILD_SLICER_EXTENSION
#       is propagated from the top-level down to this project.
if(APPLE)
  set(SlicerExecutionModel_DEFAULT_CLI_EXECUTABLE_LINK_FLAGS
    "${SlicerExecutionModel_DEFAULT_CLI_EXECUTABLE_LINK_FLAGS} -Wl,-rpath,@loader_path/../../../../../")
  set(SlicerExecutionModel_DEFAULT_CLI_LIBRARY_LINK_FLAGS
    "${SlicerExecutionModel_DEFAULT_CLI_LIBRARY_LINK_FLAGS} -Wl,-rpath,@loader_path/../../../../../")
endif()

macro(SEM_BUILD_EXECUTABLE)
  set(options )
  set( oneValueArgs NAME )
  set(multiValueArgs LIBRARIES )
  CMAKE_PARSE_ARGUMENTS(LOCAL
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
#  Useful for debugging where executable and sharedlib modules will
#  be installed.
#  message(STATUS "NAME ${LOCAL_NAME}
#    ${STATIC}
#    TARGET_LIBRARIES ${LOCAL_LIBRARIES}
#    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
#    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
#    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
#    INSTALL_RUNTIME_DESTINATION ${CLI_INSTALL_RUNTIME_DESTINATION}
#    INSTALL_LIBRARY_DESTINATION ${CLI_INSTALL_LIBRARY_DESTINATION}
#    INSTALL_ARCHIVE_DESTINATION ${CLI_INSTALL_ARCHIVE_DESTINATION}"
#  )
  SEMMacroBuildCLI(
    NAME ${LOCAL_NAME}
    ${STATIC}
    TARGET_LIBRARIES ${LOCAL_LIBRARIES}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
    INSTALL_RUNTIME_DESTINATION ${INSTALL_RUNTIME_DESTINATION}
    INSTALL_LIBRARY_DESTINATION ${INSTALL_LIBRARY_DESTINATION}
    INSTALL_ARCHIVE_DESTINATION ${INSTALL_ARCHIVE_DESTINATION}
    )
  export(TARGETS ${LOCAL_NAME} APPEND FILE ${CMAKE_BINARY_DIR}/${PRIMARY_PROJECT_NAME}-exports.cmake)
endmacro()

##dtiprocess
set( MODULE_LIBRARIES ${DTIProcess_ITK_LIBRARIES} TensorOperations DTIIO )
SEM_BUILD_EXECUTABLE( NAME dtiprocess LIBRARIES ${MODULE_LIBRARIES} )
##dtiestim
set( MODULE_LIBRARIES ${DTIProcess_ITK_LIBRARIES} cephes )
SEM_BUILD_EXECUTABLE( NAME dtiestim LIBRARIES ${MODULE_LIBRARIES} )
##fiberprocess
set( MODULE_LIBRARIES DTIIO ${DTIProcess_ITK_LIBRARIES} )
SEM_BUILD_EXECUTABLE( NAME fiberprocess LIBRARIES ${MODULE_LIBRARIES} )
##dtiaverage
set( MODULE_LIBRARIES TensorOperations DTIIO )
SEM_BUILD_EXECUTABLE( NAME dtiaverage LIBRARIES ${MODULE_LIBRARIES} )
##fiberstats
set( MODULE_LIBRARIES DTIIO )
SEM_BUILD_EXECUTABLE( NAME fiberstats LIBRARIES ${MODULE_LIBRARIES} )

#We do not build those old tools as part of the Slicer extension package. Those tools are not maintained anymore.
if( NOT DTIProcess_BUILD_SLICER_EXTENSION )
  ##scalartransform
  set( MODULE_LIBRARIES DTIIO ${DTIProcess_ITK_LIBRARIES} )
  SEM_BUILD_EXECUTABLE( NAME scalartransform LIBRARIES ${MODULE_LIBRARIES} )
  ##fibertrack
  set( MODULE_LIBRARIES DTIIO TensorOperations ${DTIProcess_ITK_LIBRARIES} )
  SEM_BUILD_EXECUTABLE( NAME fibertrack LIBRARIES ${MODULE_LIBRARIES} )
  ##maxcurvature
  set( MODULE_LIBRARIES ${DTIProcess_ITK_LIBRARIES} )
  SEM_BUILD_EXECUTABLE( NAME maxcurvature LIBRARIES ${MODULE_LIBRARIES} )
endif()

if(BUILD_dwiAtlas)
  ADD_SUBDIRECTORY(dwiAtlas)
  ADD_DEPENDENCIES(dwiAtlas dwiprocess DTIIO)
endif()


