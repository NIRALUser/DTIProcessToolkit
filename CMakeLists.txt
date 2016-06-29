project(DTIProcessToolkit)
cmake_minimum_required(VERSION 2.8.7)

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

option( DTIProcess_SUPERBUILD
        "Build DTIProcess as a Superbuild project" ON )
if( DTIProcess_SUPERBUILD )
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/DTIProcess.cmake")
endif()

