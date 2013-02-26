project(DTIProcess)
cmake_minimum_required(VERSION 2.8.7)

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

option( DTIProcess_BUILD_SLICER_EXTENSION
        "Build DTIProcess as a Slicer extension or as a Superbuild project" OFF )
if( DTIProcess_BUILD_SLICER_EXTENSION )
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/DTIProcess.cmake")
endif()

