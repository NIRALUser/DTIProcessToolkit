
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(INSTALL_RUNTIME_DESTINATION bin)
SETIFEMPTY(INSTALL_LIBRARY_DESTINATION bin)
SETIFEMPTY(INSTALL_ARCHIVE_DESTINATION lib/static)

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

find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})

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
  ITKIOTransformBase
  ITKIOTransformInsightLegacy
  ITKIOTransformMatlab
  ITKImageCompare
  ITKTestKernel
)

# IO modules: we only load the IO modules that are available with ITK. Because the toolds compiled in this project read and write DWI and DTI,
# we leave only NRRD as the compulsory module.
set( ALL_IO
  ITKIOBMP
  ITKIOBioRad
  ITKIOCSV
  ITKIODCMTK
  ITKIOGDCM
  ITKIOGE
  ITKIOGIPL
  ITKIOHDF5
  ITKIOIPL
  ITKIOJPEG
  ITKIOLSM
  ITKIOMRC
  MGHIO
  ITKIOMesh
  ITKIOMeta
  ITKIONIFTI
  ITKIOPNG
  ITKIORAW
  ITKIOSiemens
  ITKIOSpatialObjects
  ITKIOStimulate
  ITKIOTIFF
  ITKIOTransformHDF5
  ITKIOVTK
  ITKIOXML
)

list(APPEND LIST_ITK_IO_USED ITKIONRRD )
foreach( io ${ALL_IO} )
  list( FIND ITK_MODULES_ENABLED ${io} position )
  if( ${position} GREATER -1 ) #not found: ${position}==-1
    list( APPEND LIST_ITK_IO_USED ${io} )
  endif()
endforeach()

find_package(ITK COMPONENTS 
   ${LIST_ITK_COMPONENTS} ${LIST_ITK_IO_USED}
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
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/niral_utilities-install/bin/polydatamerge${fileextension} DESTINATION ${INSTALL_RUNTIME_DESTINATION} )
  endif()
  if( BUILD_PolyDataTransform )
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/niral_utilities-install/bin/polydatatransform${fileextension} DESTINATION ${INSTALL_RUNTIME_DESTINATION} )
  endif()
  if( BUILD_CropDTI )
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/niral_utilities-install/bin/CropDTI${fileextension} DESTINATION ${INSTALL_RUNTIME_DESTINATION} )
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
set( USE_SYSTEM_ITK ON)
set( USE_SYSTEM_VTK ON)
set( USE_SYSTEM_SlicerExecutionModel ON)

set(extProjName ${PRIMARY_PROJECT_NAME})
set(proj        ${PRIMARY_PROJECT_NAME})
List( LENGTH ${PRIMARY_PROJECT_NAME}_DEPENDENCIES dependencies_size )
if( dependencies_size GREATER 0 )
  SlicerMacroCheckExternalProjectDependency(${PRIMARY_PROJECT_NAME})
endif()

IF(BUILD_TESTING)
  include(CTest)
  ADD_SUBDIRECTORY(Testing)
ENDIF(BUILD_TESTING)


if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/DTIProcess)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

configure_file(CMake/DTIProcessConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/DTIProcessConfig.cmake" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/DTIProcessConfig.cmake"  
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)