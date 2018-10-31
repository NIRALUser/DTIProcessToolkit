


if( DTIProcess_BUILD_SLICER_EXTENSION )
  unsetForSlicer( NAMES QT_QMAKE_EXECUTABLE SlicerExecutionModel_DIR ITK_DIR VTK_DIR CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS CMAKE_C_FLAGS ITK_LIBRARIES )
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  resetForSlicer( NAMES CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS CMAKE_C_FLAGS )
  
  SET(INSTALL_RUNTIME_DESTINATION ${Slicer_INSTALL_CLIMODULES_BIN_DIR})
  SET(INSTALL_LIBRARY_DESTINATION ${Slicer_INSTALL_CLIMODULES_LIB_DIR})
  SET(INSTALL_ARCHIVE_DESTINATION ${Slicer_INSTALL_CLIMODULES_LIB_DIR})
  
endif()

SETIFEMPTY(INSTALL_RUNTIME_DESTINATION bin)
SETIFEMPTY(INSTALL_LIBRARY_DESTINATION bin)
SETIFEMPTY(INSTALL_ARCHIVE_DESTINATION lib/static)
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

option(BUILD_TESTING "Build the testing tree" ON)
option(EXECUTABLES_ONLY "Build only executables (CLI)" OFF)

if( ${EXECUTABLES_ONLY} )
  set( STATIC "EXECUTABLE_ONLY" )
  set( STATIC_LIB "STATIC" )
else()
  set( STATIC_LIB "SHARED" )
endif()


find_package(niral_utilities)
if(niral_utilities_FOUND)

  foreach(niral_utilities_lib ${niral_utilities_LIBRARIES})

    get_target_property(niral_utilities_location ${niral_utilities_lib} LOCATION_RELEASE)
    if(NOT EXISTS ${niral_utilities_location})
      message(STATUS "skipping niral_utilities_lib install rule: [${niral_utilities_location}] does not exist")
      continue()
    endif()
    
    if(EXISTS "${niral_utilities_location}.xml")
      install(PROGRAMS ${niral_utilities_location} 
        DESTINATION ${INSTALL_RUNTIME_DESTINATION}
        COMPONENT RUNTIME)

      install(FILES ${niral_utilities_location}.xml
        DESTINATION ${INSTALL_RUNTIME_DESTINATION}
        COMPONENT RUNTIME)
    else()
      install(PROGRAMS ${niral_utilities_location} 
        DESTINATION ${INSTALL_RUNTIME_DESTINATION}/../ExternalBin
        COMPONENT RUNTIME)      
    endif()
  endforeach()
  
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


find_package(VTK COMPONENTS
  vtkIOLegacy
  vtkIOXML
  vtkCommonDataModel
  REQUIRED)
include(${VTK_USE_FILE})



INCLUDE_DIRECTORIES(
${PROJECT_SOURCE_DIR}/Library
${PROJECT_SOURCE_DIR}/PrivateLibrary
${PROJECT_SOURCE_DIR}
)


## Replace bessel(FORTRAN) with cephes(C)
SET(BESSEL_LIB cephes)
ADD_SUBDIRECTORY(cephes)

ADD_SUBDIRECTORY(PrivateLibrary)
ADD_SUBDIRECTORY(Applications)

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

if( DTIProcess_BUILD_SLICER_EXTENSION )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ImportDTIProcessExtensionExecutables.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ImportDTIProcessExtensionExecutables.cmake)
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()

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

if( DTIAtlasFiberAnalyzer_BUILD_SLICER_EXTENSION )
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()