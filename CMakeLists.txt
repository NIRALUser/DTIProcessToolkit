
cmake_minimum_required(VERSION 2.8.7)

option( DTIProcess_SUPERBUILD
        "Build DTIProcess as a Superbuild project" ON )
if( DTIProcess_SUPERBUILD )
  project(DTIProcessToolkit)
  include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
else()
  project(DTIProcess)
  include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)
  include("${CMAKE_CURRENT_SOURCE_DIR}/DTIProcess.cmake")
endif()

