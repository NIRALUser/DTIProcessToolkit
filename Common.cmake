include(CMakeDependentOption)

enable_language(C)
enable_language(CXX)

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
if( DTIProcess_BUILD_SLICER_EXTENSION )
  find_package(Subversion REQUIRED)
endif()

find_package(Git)
if(NOT GIT_FOUND)
  message(WARNING "Git may be needed to download external dependencies: Install Git and try to re-configure")
endif()

option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
else(NOT USE_GIT_PROTOCOL)
  set(git_protocol "git")
endif()

#-----------------------------------------------------------------------------
# Platform check
#-----------------------------------------------------------------------------

set(PLATFORM_CHECK true)

if(PLATFORM_CHECK)
  # See CMake/Modules/Platform/Darwin.cmake)
  #   6.x == Mac OSX 10.2 (Jaguar)
  #   7.x == Mac OSX 10.3 (Panther)
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  if (DARWIN_MAJOR_VERSION LESS "9")
    message(FATAL_ERROR "Only Mac OSX >= 10.5 are supported !")
  endif()
endif()

#-----------------------------------------------------------------------------
set(EXTENSION_NAME DTIProcess)
set(EXTENSION_HOMEPAGE "http://www.nitrc.org/projects/dtiprocess")
set(EXTENSION_CATEGORY "Diffusion")
set(EXTENSION_CONTRIBUTORS "Francois Budin (UNC)")
set(EXTENSION_DESCRIPTION "This extension provides the tool DTIProcess integrated in Slicer")
set(EXTENSION_ICONURL "http://www.nitrc.org/project/screenshot.php?group_id=312&screenshot_id=771")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/thumb/b/b8/DTIEstim-B0-crop.png/193px-DTIEstim-B0-crop.png http://www.slicer.org/slicerWiki/images/thumb/9/90/FiberTrack-fibers.png/138px-FiberTrack-fibers.png")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated list or 'NA' if any
set(EXTENSION_BUILD_SUBDIRECTORY DTIProcess-build)


option( BUILD_TESTING   "Build the testing tree" ON )


# With CMake 2.8.9 or later, the UPDATE_COMMAND is required for updates to occur.
# For earlier versions, we nullify the update state to prevent updates and
# undesirable rebuild.
if(CMAKE_VERSION VERSION_LESS 2.8.9)
  set(cmakeversion_external_update UPDATE_COMMAND "")
else()
  set(cmakeversion_external_update LOG_UPDATE 1)
endif()

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(SlicerExtensionsConfigureMacros)
#-----------------------------------------------------------------------------
# CMake Function(s) and Macro(s)
#-----------------------------------------------------------------------------
if(CMAKE_VERSION VERSION_LESS 2.8.3)
  include(Pre283CMakeParseArguments)
else()
  include(CMakeParseArguments)
endif()

#-----------------------------------------------------------------------------
if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

#-------------------------------------------------------------------------
SETIFEMPTY(CLI_INSTALL_LIBRARY_DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SETIFEMPTY(CLI_INSTALL_ARCHIVE_DESTINATION ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
SETIFEMPTY(CLI_INSTALL_RUNTIME_DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


#-------------------------------------------------------------------------
# Augment compiler flags
#-------------------------------------------------------------------------
include(ITKSetStandardCompilerFlags)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS}" )
else() # Release, or anything else
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS}" )
endif()

#-----------------------------------------------------------------------------
# Add needed flag for gnu on linux like enviroments to build static common libs
# suitable for linking with shared object libs.
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
  if(NOT "${CMAKE_C_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  endif()
endif()
