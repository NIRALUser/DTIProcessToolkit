PROJECT(DTIProcess)
cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)
#INCLUDE(CPack)

ENABLE_TESTING()
INCLUDE(Dart)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DTIProcess_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DTIProcess_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DTIProcess_BINARY_DIR}/bin)

IF( BUILD_AGAINST_SLICER3 )
  find_package(Slicer3 REQUIRED)
  include(${Slicer3_USE_FILE})

  slicer3_set_default_install_prefix_for_external_projects()
ENDIF( BUILD_AGAINST_SLICER3 )


FIND_PACKAGE( ITK REQUIRED )
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE( VTK REQUIRED )
INCLUDE(${VTK_USE_FILE})

##  In many cases sub-projects depending on SlicerExectuion Model
##  that can be built stand alone are combined in larger packages
##  This logic will include SlicerExectionModel only if it
##  has not already been built.
FIND_PACKAGE(GenerateCLP NO_MODULE QUIET)
IF( NOT GenerateCLP_DIR )
      ADD_SUBDIRECTORY(SlicerExecutionModel)
ENDIF( NOT GenerateCLP_DIR )
FIND_PACKAGE(GenerateCLP NO_MODULE REQUIRED)
IF(GenerateCLP_DIR)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_DIR)
  MESSAGE(FATAL_ERROR "Can't build without GenerateCLP. Please set GenerateCLP_DIR")
ENDIF(GenerateCLP_DIR)


INCLUDE_DIRECTORIES(
${DTIProcess_SOURCE_DIR}/Library
${DTIProcess_SOURCE_DIR}/PrivateLibrary
${DTIProcess_SOURCE_DIR}
)

## Replace bessel(FORTRAN) with cephes(C)
SET(BESSEL_LIB cephes)
ADD_SUBDIRECTORY(cephes)

ADD_SUBDIRECTORY(Library)
ADD_SUBDIRECTORY(PrivateLibrary)
ADD_SUBDIRECTORY(Applications)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(Testing)
ENDIF(BUILD_TESTING)
