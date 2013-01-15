cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  set(EXTENSION_NAME DTIProcess)
  set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.2/Extensions/DTIProcess")
  set(EXTENSION_CATEGORY "Diffusion")
  set(EXTENSION_CONTRIBUTORS "Francois Budin (UNC)")
  set(EXTENSION_DESCRIPTION "This extension provides the tool DTIProcess integrated in Slicer")
  set(EXTENSION_ICONURL "http://www.nitrc.org/project/screenshot.php?group_id=312&screenshot_id=575")
  set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/thumb/b/b8/DTIEstim-B0-crop.png/193px-DTIEstim-B0-crop.png http://www.slicer.org/slicerWiki/images/thumb/9/90/FiberTrack-fibers.png/138px-FiberTrack-fibers.png")
  set(EXTENSION_STATUS "Beta")
  set(EXTENSION_DEPENDS "NA") # Specified as a space separated list or 'NA' if any
  set(EXTENSION_BUILD_SUBDIRECTORY dtiprocess-build)
endif()

#-----------------------------------------------------------------------------
set(MODULE_NAME DTIProcessExtension)

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})




include(ExternalProject)
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()



option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
else(NOT USE_GIT_PROTOCOL)
  set(git_protocol "git")
endif()

set(proj Insight)
set(${proj}_REPOSITORY ${git_protocol}://itk.org/ITK.git CACHE STRING "" FORCE)
ExternalProject_Add(${proj}
  GIT_TAG v4.3.0
  GIT_REPOSITORY ${${proj}_REPOSITORY}
  UPDATE_COMMAND ""
  SOURCE_DIR ${proj}
  BINARY_DIR ${proj}-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DITKV3_COMPATIBILITY:BOOL=ON
    -DBUILD_TESTING:BOOL=OFF
    -DBUILD_EXAMPLES:BOOL=OFF
  INSTALL_COMMAND ""
  )
set(ITK4_DIR ${CMAKE_CURRENT_BINARY_DIR}/Insight-build)
set(ITK4_DEPEND ${proj}) ## Set the internal dependancy for ITK


set(proj SlicerExecutionModel)
  ExternalProject_Add(${proj}
  GIT_REPOSITORY ${git_protocol}://github.com/Slicer/SlicerExecutionModel.git
  GIT_TAG "5ac91a89bba50db8b4f5e32cadbcb4baadb1ffc0"
  SOURCE_DIR ${proj}
  BINARY_DIR ${proj}-build
  DEPENDS ${ITK4_DEPEND}
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DITK_DIR:PATH=${ITK4_DIR}
    INSTALL_COMMAND ""
  )

set(SlicerExecutionModel4_DIR ${CMAKE_CURRENT_BINARY_DIR}/SlicerExecutionModel-build)
set(SlicerExecutionModel4_DEPEND ${proj}) ## Set the internal dependancy for ITK

  set(proj dtiprocess)
  ExternalProject_Add(${proj}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dtiprocess
    BINARY_DIR ${proj}-build
    DEPENDS  ${ITK4_DEPEND} ${SlicerExecutionModel4_DEPEND} ${VTK_DEPEND}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${LOCAL_CMAKE_BUILD_OPTIONS}
      -DBUILD_TESTING:BOOL=ON
      -DITK_DIR:PATH=${ITK4_DIR}
      -DVTK_DIR:PATH=${VTK_DIR}
      -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel4_DIR}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
      -DCLI_RUNTIME_OUTPUT_DIRECTORY:PATH=${SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY}
      -DCLI_INSTALL_RUNTIME_DESTINATION:PATH=${SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION}
      -DMIDAS_PACKAGE_EMAIL:STRING=${MIDAS_PACKAGE_EMAIL}
      -DMIDAS_PACKAGE_API_KEY:STRING=${MIDAS_PACKAGE_API_KEY}
      -DEXTENSION_NAME:STRING=${EXTENSION_NAME}
      -DEXTENSION_SUPERBUILD_BINARY_DIR:PATH=${${EXTENSION_NAME}_BINARY_DIR}
      # Slicer
      -DSlicer_DIR:PATH=${Slicer_DIR}
      -DSlicer_USE_FILE:FILEPATH=${Slicer_USE_FILE}
    INSTALL_COMMAND ""
  )

#-----------------------------------------------------------------------------
