#-----------------------------------------------------------------------------
set(MODULE_NAME ${EXTENSION_NAME}) # Do not use 'project()'
set(MODULE_TITLE ${MODULE_NAME})

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)
unset( USE_SYSTEM_ITK CACHE )
unset( USE_SYSTEM_VTK CACHE )
unset( USE_SYSTEM_SlicerExecutionModel CACHE )

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )


## A simple macro to set variables ONLY if it has not been set
## This is needed when stand-alone packages are combined into
## a larger package, and the desired behavior is that all the
## binary results end up in the combined directory.
if(NOT SETIFEMPTY)
macro(SETIFEMPTY)
  set(KEY ${ARGV0})
  set(VALUE ${ARGV1})
  if(NOT ${KEY})
    set(${KEY} ${VALUE})
  endif(NOT ${KEY})
endmacro(SETIFEMPTY KEY VALUE)
endif(NOT SETIFEMPTY)
###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(INSTALL_RUNTIME_DESTINATION bin)
SETIFEMPTY(INSTALL_LIBRARY_DESTINATION bin)
SETIFEMPTY(INSTALL_ARCHIVE_DESTINATION lib/static)

option(BUILD_dwiAtlas "Build dwiAtlas or not.  Requires boost." OFF)
option(BUILD_TESTING "Build the testing tree" ON)
option(EXECUTABLES_ONLY "Build only executables (CLI)" OFF)

if( ${EXECUTABLES_ONLY} )
  set( STATIC "EXECUTABLE_ONLY" )
  set( STATIC_LIB "STATIC" )
else()
  set( STATIC_LIB "SHARED" )
endif()

##  In many cases sub-projects depending on SlicerExectuion Model
##  that can be built stand alone are combined in larger packages.
##  This logic will include SlicerExectionModel only if it
##  has not already been included by a previous package.

if( DTIProcess_BUILD_SLICER_EXTENSION )
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()


find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)
if( NOT DEFINED ITKV3_COMPATIBILITY OR NOT ${ITKV3_COMPATIBILITY}  )
  message( WARNING "Choose ITKv4 compiled with ITKV3_COMPATIBILITY set to ON (or GenerateCLP compiled against such an ITK version). If not, you may have compilation errors" )
endif()



if(NOT VTK_FOUND)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif(NOT VTK_FOUND)


INCLUDE_DIRECTORIES(
${DTIProcess_SOURCE_DIR}/Library
${DTIProcess_SOURCE_DIR}/PrivateLibrary
${DTIProcess_SOURCE_DIR}
)


## Replace bessel(FORTRAN) with cephes(C)
SET(BESSEL_LIB cephes)
ADD_SUBDIRECTORY(cephes)

ADD_SUBDIRECTORY(PrivateLibrary)
ADD_SUBDIRECTORY(Applications)




IF(BUILD_TESTING)
  include(CTest)
  ADD_SUBDIRECTORY(Testing)
ENDIF(BUILD_TESTING)


if( DTIProcess_BUILD_SLICER_EXTENSION )
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
else()
  if( NOT WIN32 )
    set( Tools  dtiaverage dtiestim dtiprocess fiberprocess fiberstats fibertrack maxcurvature scalartransform )
    foreach( tool ${Tools} )
      install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${tool} DESTINATION ${INSTALL_RUNTIME_DESTINATION})
    endforeach()
  else()
    message( WARNING "No install on Windows" )
  endif()
endif()
