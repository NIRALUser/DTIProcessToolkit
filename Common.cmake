#-----------------------------------------------------------------------------
set(PRIMARY_PROJECT_NAME DTIProcess)

#-----------------------------------------------------------------------------
# Standalone vs Slicer extension option
#-----------------------------------------------------------------------------

# This option should be named after the project name, it corresponds to the
# option set to ON when the project is build by the Slicer Extension build
# system.

set(_default OFF)
set(_reason "${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION is ON")
if(NOT DEFINED ${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION AND DEFINED Slicer_DIR)
  set(_default ON)
  set(_reason "Slicer_DIR is SET")
endif()

option(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION "Build as a Slicer Extension" ${_default})

set(_msg "Checking if building as a Slicer extension")
message(STATUS ${_msg})
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  message(STATUS "${_msg} - yes (${_reason})")
else()
  message(STATUS "${_msg} - no (${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION is OFF)")
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
set(verbose FALSE)
#-----------------------------------------------------------------------------
if(WIN32)
  set(fileextension .exe)
endif()
#-----------------------------------------------------------------------------

include(CMakeDependentOption)

enable_language(C)
enable_language(CXX)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------
include(ExternalProject)
include(SlicerMacroEmptyExternalProject)
include(SlicerMacroCheckExternalProjectDependency)

if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
find_package(Subversion REQUIRED)

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
set(EXTENSION_BUILD_SUBDIRECTORY .)
set(EXTENSION_README_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
set(EXTENSION_LICENSE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/License.txt)

option( BUILD_TESTING   "Build the testing tree" ON )


#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(SlicerExtensionsConfigureMacros)
#-----------------------------------------------------------------------------
# CMake Function(s) and Macro(s)
#-----------------------------------------------------------------------------
include(CMakeParseArguments)

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

option(BUILD_dwiAtlas "Build dwiAtlas or not.  Requires boost." OFF)
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

# Since C++11 support should NOT be enabled when building against Slicer 4.8, we
# use the following workaround to find out which version of Slicer is used.
# Indeed, simply calling "find_package(Slicer ...)" here would have undesirable side
# effects. For more details, see use of "unsetForSlicer()" in SuperBuild.cmake
if(NOT DEFINED HAS_STD_CPP11_FLAG)
  set(_msg "Forcing HAS_STD_CPP11_FLAG to OFF")
  message(STATUS "${_msg}")
  if(DEFINED Slicer_DIR)
    set(_config "${Slicer_DIR}/SlicerConfig.cmake")
    if(NOT EXISTS ${_config})
      message(FATAL_ERROR "Slicer_DIR [${Slicer_DIR}] is specified but doesn't contain SlicerConfig.cmake")
    endif()
    file(READ "${Slicer_DIR}/SlicerConfig.cmake" _config_content)
    string(REGEX REPLACE "Slicer_REQUIRED_QT_VERSION \"([0-9\\.]+)\"" "\\1" _match ${_config_content})
    if(_match STREQUAL "")
      message(FATAL_ERROR "Failed to extract Slicer_REQUIRED_QT_VERSION")
    endif()
    set(Slicer_REQUIRED_QT_VERSION "${CMAKE_MATCH_1}")
    if(Slicer_REQUIRED_QT_VERSION VERSION_LESS 5)
      set(HAS_STD_CPP11_FLAG 0)
      set(_status "yes (Slicer built using Qt4 requires C++98)")
    else()
      set(_status "no (Slicer built using Qt5 supports C++11)")
    endif()
  else()
    set(_status "no (Standalone build of DTIProcess supports C++11)")
  endif()
  message(STATUS "${_msg} - ${_status}")
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++11 HAS_STD_CPP11_FLAG)
if(HAS_STD_CPP11_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()

#-----------------------------------------------------------------------------
# Define Superbuild global variables
#-----------------------------------------------------------------------------
set(EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Select where external packages will be downloaded" )
set(EXTERNAL_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Select where external packages will be compiled and installed" )
# This variable will contain the list of CMake variable specific to each external project
# that should passed to ${CMAKE_PROJECT_NAME}.
# The item of this list should have the following form: <EP_VAR>:<TYPE>
# where '<EP_VAR>' is an external project variable and TYPE is either BOOL, STRING, PATH or FILEPATH.
# TODO Variable appended to this list will be automatically exported in ${PRIMARY_PROJECT_NAME}Config.cmake,
# prefix '${PRIMARY_PROJECT_NAME}_' will be prepended if it applies.
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS)

# The macro '_expand_external_project_vars' can be used to expand the list of <EP_VAR>.
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS) # List of CMake args to configure BRAINS
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES) # List of CMake variable names

# Convenient macro allowing to expand the list of EP_VAR listed in ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
# The expanded arguments will be appended to the list ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS
# Similarly the name of the EP_VARs will be appended to the list ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES.
macro(_expand_external_project_vars)
  set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS "")
  set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES "")
  foreach(arg ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS})
    string(REPLACE ":" ";" varname_and_vartype ${arg})
    set(target_info_list ${target_info_list})
    list(GET varname_and_vartype 0 _varname)
    list(GET varname_and_vartype 1 _vartype)
    list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS -D${_varname}:${_vartype}=${${_varname}})
    list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES ${_varname})
  endforeach()
endmacro()

#-----------------------------------------------------------------------------
# Common external projects CMake variables
#-----------------------------------------------------------------------------
list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
  MAKECOMMAND:STRING
  CMAKE_SKIP_RPATH:BOOL
  CMAKE_MODULE_PATH:PATH
  CMAKE_BUILD_TYPE:STRING
  BUILD_SHARED_LIBS:BOOL
  CMAKE_CXX_COMPILER:PATH
  CMAKE_CXX_FLAGS_RELEASE:STRING
  CMAKE_CXX_FLAGS_DEBUG:STRING
  CMAKE_CXX_FLAGS:STRING
  CMAKE_C_COMPILER:PATH
  CMAKE_C_FLAGS_RELEASE:STRING
  CMAKE_C_FLAGS_DEBUG:STRING
  CMAKE_C_FLAGS:STRING
  #CMAKE_CXX_STANDARD:STRING
  CMAKE_SHARED_LINKER_FLAGS:STRING
  CMAKE_EXE_LINKER_FLAGS:STRING
  CMAKE_MODULE_LINKER_FLAGS:STRING
  CMAKE_GENERATOR:STRING
  CMAKE_EXTRA_GENERATOR:STRING
  CMAKE_INSTALL_PREFIX:PATH
  CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
  CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
  CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
  CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
  CTEST_NEW_FORMAT:BOOL
  MEMORYCHECK_COMMAND_OPTIONS:STRING
  MEMORYCHECK_COMMAND:PATH
  CMAKE_SHARED_LINKER_FLAGS:STRING
  CMAKE_EXE_LINKER_FLAGS:STRING
  CMAKE_MODULE_LINKER_FLAGS:STRING
  SITE:STRING
  BUILDNAME:STRING
  Subversion_SVN_EXECUTABLE:FILEPATH
  GIT_EXECUTABLE:FILEPATH
  USE_GIT_PROTOCOL:BOOL
  Slicer_DIR:PATH
  VTK_VERSION_MAJOR:STRING
  )

option(BUILD_PolyDataTransform "Build PolyDataTransform" ON)
option(BUILD_PolyDataMerge "Build PolyDataMerge" ON)
option(BUILD_CropDTI "Build CropDTI" ON)

