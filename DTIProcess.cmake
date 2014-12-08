#-----------------------------------------------------------------------------
set(MODULE_NAME ${EXTENSION_NAME}) # Do not use 'project()'
set(MODULE_TITLE ${MODULE_NAME})

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)
#unset( USE_SYSTEM_ITK CACHE )
#unset( USE_SYSTEM_VTK CACHE )
#unset( USE_SYSTEM_SlicerExecutionModel CACHE )


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

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

set( LIST_ITK_COMPONENTS
  ITKDiffusionTensorImage
  ITKRegistrationCommon
  ITKOptimizersv4
  ITKConnectedComponents
  ITKMathematicalMorphology
  ITKBinaryMathematicalMorphology
  ITKRegionGrowing
  ITKMetricsv4
  ITKRegistrationMethodsv4
  ITKDistanceMap
  ITKVTK
  ITKTransform
  ITKIOImageBase
  ITKIONRRD
  ITKImageCompare
  ITKIOBMP
  ITKIOBioRad
  ITKIOCSV
  ITKIODCMTK
  ITKIOGDCM
  ITKIOGE
  ITKIOGIPL
  ITKIOHDF5
  ITKIOIPL
  ITKIOImageBase
  ITKIOJPEG
  ITKIOLSM
  ITKIOMRC
  ITKIOMesh
  ITKIOMeta
  ITKIONIFTI
  ITKIONRRD
  ITKIOPNG
  ITKIORAW
  ITKIOSiemens
  ITKIOSpatialObjects
  ITKIOStimulate
  ITKIOTIFF
  ITKIOTransformBase
  ITKIOTransformHDF5
  ITKIOTransformInsightLegacy
  ITKIOTransformMatlab
  ITKIOVTK
  ITKIOXML
  ITKTestKernel
)

if( NOT DTIProcess_BUILD_SLICER_EXTENSION )
  list( APPEND LIST_ITK_COMPONENTS MGHIO )
endif()

find_package(ITK COMPONENTS 
   ${LIST_ITK_COMPONENTS}
   REQUIRED
   )
include(${ITK_USE_FILE})
# have to save library list because it gets clobbered by GenerateCLP, which
# calls find_package(ITK) every time it's invoked.
set(DTIProcess_ITK_LIBRARIES ${ITK_LIBRARIES})



if(NOT VTK_FOUND)
    find_package(VTK COMPONENTS
      vtkIOLegacy
      vtkIOXML
      REQUIRED)
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


option(BUILD_PolyDataTransform "Build PolyDataTransform" ON)
option(BUILD_PolyDataMerge "Build PolyDataMerge" ON)
option(BUILD_CropDTI "Build CropDTI" ON)
if( BUILD_PolyDataTransform OR BUILD_PolyDataMerge OR BUILD_CropDTI )
  set(${PRIMARY_PROJECT_NAME}_DEPENDENCIES niral_utilities)
  if( BUILD_PolyDataMerge )
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/niral_utilities-install/bin/polydatamerge DESTINATION ${INSTALL_RUNTIME_DESTINATION} )
  endif()
  if( BUILD_PolyDataTransform )
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/niral_utilities-install/bin/polydatatransform DESTINATION ${INSTALL_RUNTIME_DESTINATION} )
  endif()
  if( BUILD_CropDTI )
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/niral_utilities-install/bin/CropDTI DESTINATION ${INSTALL_RUNTIME_DESTINATION} )
  endif()
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
set(extProjName ${PRIMARY_PROJECT_NAME})
set(proj        ${PRIMARY_PROJECT_NAME})
List( LENGTH ${PRIMARY_PROJECT_NAME}_DEPENDENCIES dependencies_size )
if( dependencies_size GREATER 0 )
  SlicerMacroCheckExternalProjectDependency(${PRIMARY_PROJECT_NAME})
endif()


if( DTIProcess_BUILD_SLICER_EXTENSION )
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

IF(BUILD_TESTING)
  include(CTest)
  ADD_SUBDIRECTORY(Testing)
ENDIF(BUILD_TESTING)


if( DTIProcess_BUILD_SLICER_EXTENSION )
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()
