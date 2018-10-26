set(PRIMARY_PROJECT_NAME DTIProcess)

#-----------------------------------------------------------------------------
# Superbuild option(s)
#-----------------------------------------------------------------------------

option(EXECUTABLES_ONLY "Build the tools and the tools' libraries statically" ON)
#-----------------------------------------------------------------------------
SETIFEMPTY(INSTALL_RUNTIME_DESTINATION bin)
SETIFEMPTY(INSTALL_LIBRARY_DESTINATION lib)
set( DTIProcess_INSTALL_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/DTIProcess-install )
#-----------------------------------------------------------------------------
# Extension configuration
#-----------------------------------------------------------------------------
########Depending if it is an extension or a Superbuild
if( DTIProcess_BUILD_SLICER_EXTENSION )
  unsetForSlicer(NAMES CMAKE_MODULE_PATH CMAKE_C_COMPILER CMAKE_CXX_COMPILER DCMTK_DIR ITK_DIR SlicerExecutionModel_DIR VTK_DIR QT_QMAKE_EXECUTABLE ITK_VERSION_MAJOR CMAKE_CXX_FLAGS CMAKE_C_FLAGS Teem_DIR)
  find_package(Slicer REQUIRED)
  set( Slicer_USE_PYTHONQT FALSE )
  set( USE_SYSTEM_ITK ON CACHE BOOL "Build using an externally defined version of ITK" FORCE )
  set( USE_SYSTEM_VTK ON CACHE BOOL "Build using an externally defined version of VTK" FORCE )
  #VTK_VERSION_MAJOR is define but not a CACHE variable
  set( VTK_VERSION_MAJOR ${VTK_VERSION_MAJOR} CACHE STRING "Choose the expected VTK major version to build Slicer (5 or 6).")
  set( USE_SYSTEM_SlicerExecutionModel ON CACHE BOOL "Build using an externally defined version of SlicerExecutionModel" FORCE )
  set( BUILD_CropDTI OFF CACHE BOOL "Build CropDTI" FORCE) # CropDTI is not a CLI
  set(COMPILE_IMAGEMATH ON)

  SET(INSTALL_RUNTIME_DESTINATION ${Slicer_INSTALL_CLIMODULES_BIN_DIR})
  SET(INSTALL_LIBRARY_DESTINATION ${Slicer_INSTALL_CLIMODULES_LIB_DIR})
  SET(INSTALL_ARCHIVE_DESTINATION ${Slicer_INSTALL_CLIMODULES_LIB_DIR})
endif()

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)

#------------------------------------------------------------------------------
# ${PRIMARY_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------
#Compilation options for ITK if USE_SYSTEM_ITK is set to OFF (Superbuild)
set(ITK_EXTERNAL_NAME ITKv4)
#set(${PRIMARY_PROJECT_NAME}_BUILD_DICOM_SUPPORT ON)
set( USE_ITK_Module_MGHIO ON )
if( UNIX )
  set( ${PROJECT_NAME}_BUILD_TIFF_SUPPORT ON )
  set( ${PROJECT_NAME}_BUILD_JPEG_SUPPORT ON )
endif()
set( ${PRIMARY_PROJECT_NAME}_BUILD_ZLIB_SUPPORT ON )
set( ${PRIMARY_PROJECT_NAME}_BUILD_FFTW_SUPPORT ON )

set(${PRIMARY_PROJECT_NAME}_DEPENDENCIES ${ITK_EXTERNAL_NAME} SlicerExecutionModel VTK )
if( BUILD_dwiAtlas )
  list(APPEND ${PRIMARY_PROJECT_NAME}_DEPENDENCIES Boost)
endif()


option(BUILD_PolyDataTransform "Build PolyDataTransform from niral_utilities" ON)
option(BUILD_PolyDataMerge "Build PolyDataMerge from niral_utilities" ON)
option(BUILD_CropDTI "Build CropDTI from niral_utilities" ON)

if( BUILD_PolyDataTransform OR BUILD_PolyDataMerge OR BUILD_CropDTI )
  set(${PRIMARY_PROJECT_NAME}_DEPENDENCIES  ${${PRIMARY_PROJECT_NAME}_DEPENDENCIES} niral_utilities)
endif()

CMAKE_DEPENDENT_OPTION( ITK_LEGACY_REMOVE "Remove ITK legacy" ON "NOT BUILD_dwiAtlas" OFF)
CMAKE_DEPENDENT_OPTION( ITKV3_COMPATIBILITY "Build ITKv4 with ITKv3 compatibility" OFF "NOT BUILD_dwiAtlas" ON)

if(${PRIMARY_PROJECT_NAME}_USE_QT)
  list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
    ${PRIMARY_PROJECT_NAME}_USE_QT:BOOL
    QT_QMAKE_EXECUTABLE:PATH
    QT_MOC_EXECUTABLE:PATH
    QT_UIC_EXECUTABLE:PATH
    )
endif()

_expand_external_project_vars()
set(COMMON_EXTERNAL_PROJECT_ARGS ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})
set(extProjName ${PRIMARY_PROJECT_NAME})
set(proj        ${PRIMARY_PROJECT_NAME})

List( LENGTH ${PRIMARY_PROJECT_NAME}_DEPENDENCIES dependencies_size )
if( dependencies_size GREATER 0 )
  SlicerMacroCheckExternalProjectDependency(${PRIMARY_PROJECT_NAME})
endif()
#-----------------------------------------------------------------------------
# Set CMake OSX variable to pass down the external project
#-----------------------------------------------------------------------------
set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
if(APPLE)
  list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
endif()


#-----------------------------------------------------------------------------
# Add external project CMake args
#-----------------------------------------------------------------------------
list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
  BUILD_TESTING:BOOL
  ITK_VERSION_MAJOR:STRING
  ITK_DIR:PATH
  VTK_DIR:PATH
  GenerateCLP_DIR:PATH
  SlicerExecutionModel_DIR:PATH
  )

_expand_external_project_vars()
set(COMMON_EXTERNAL_PROJECT_ARGS ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})


if(verbose)
  message("Inner external project args:")
  foreach(arg ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})
    message("  ${arg}")
  endforeach()
endif()

string(REPLACE ";" "^" ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES "${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES}")

if(verbose)
  message("Inner external project argnames:")
  foreach(argname ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES})
    message("  ${argname}")
  endforeach()
endif()

#------------------------------------------------------------------------------
# Configure and build
#------------------------------------------------------------------------------

SET(CMAKE_ARGS 
      -DDTIProcess_SUPERBUILD:BOOL=OFF
      -DDTIProcess_BUILD_SLICER_EXTENSION:BOOL=${DTIProcess_BUILD_SLICER_EXTENSION}
      -DDTIProcess_EXTENSION:BOOL=ON #install the tests if it is built as an extension
      -DINSTALL_RUNTIME_DESTINATION:PATH=${INSTALL_RUNTIME_DESTINATION}
      -DINSTALL_LIBRARY_DESTINATION:PATH=${INSTALL_LIBRARY_DESTINATION}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DEXECUTABLES_ONLY:BOOL=${EXECUTABLES_ONLY}
      -DUSE_SYSTEM_ITK:BOOL=TRUE
      -DUSE_SYSTEM_VTK:BOOL=TRUE
      -DUSE_SYSTEM_SlicerExecutionModel:BOOL=TRUE
      -DBUILD_dwiAtlas:BOOL=${BUILD_dwiAtlas}
      -DCMAKE_INSTALL_PREFIX:PATH=${DTIProcess_INSTALL_DIRECTORY}
      -DCMAKE_PREFIX:PATH=${niral_utilities_DIR}
  )

set(proj DTIProcess-inner)
ExternalProject_Add(${proj}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    BINARY_DIR ${proj}-build
    DEPENDS ${${PRIMARY_PROJECT_NAME}_DEPENDENCIES}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${CMAKE_ARGS}
  )

if( DTIProcess_BUILD_SLICER_EXTENSION )
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ImportDTIProcessExtensionExecutables.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ImportDTIProcessExtensionExecutables.cmake)
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR}/DTIProcess-inner-build;${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()

#-----------------------------------------------------------------------------

