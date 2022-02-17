
set(proj ITK)
# Set dependency list
set(${proj}_DEPENDENCIES "zlib")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(ITK_DIR CACHE)
  find_package(ITK 5.3 COMPONENTS ${${CMAKE_PROJECT_NAME}_ITK_COMPONENTS} REQUIRED NO_MODULE)
endif()

# Sanity checks
if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED ITK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/Slicer/ITK.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "be81e6223240508642b963511e6441203df6375e" # Version of ITK to be used in Slicer post 4.13 release
    QUIET
    )

  if(NOT ${CMAKE_PROJECT_NAME}ITKV3_COMPATIBILITY AND CMAKE_CL_64)
    # enables using long long type for indexes and size on platforms
    # where long is only 32-bits (msvc)
    set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DITK_USE_64BITS_IDS:BOOL=ON
      )
  endif()

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    #GIT_SHALLOW TRUE
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    CMAKE_CACHE_ARGS
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      -DITK_INSTALL_ARCHIVE_DIR:PATH=${Slicer_INSTALL_LIB_DIR}
      -DITK_INSTALL_LIBRARY_DIR:PATH=${Slicer_INSTALL_LIB_DIR}
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DITK_LEGACY_REMOVE:BOOL=OFF
      -DITK_FUTURE_LEGACY_REMOVE:BOOL=OFF
      #-DITK_LEGACY_REMOVE:BOOL=ON
      #-DITK_FUTURE_LEGACY_REMOVE:BOOL=ON
      -DITKV3_COMPATIBILITY:BOOL=OFF
      -DITKVTK:BOOL=ON
      -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF
      -DModule_ITKBinaryMathematicalMorphology:BOOL=ON
      -DModule_ITKCommon:BOOL=ON
      -DModule_ITKConnectedComponents:BOOL=ON
      -DModule_ITKDeprecated:BOOL=ON # Needed for WarpImageFilter and VectorInterpolators
      -DModule_ITKDiffusionTensorImage:BOOL=ON # For niral_utilities
      -DModule_ITKDistanceMap:BOOL=ON
      -DModule_ITKIOBMP:BOOL=ON
      -DModule_ITKIOBioRad:BOOL=ON
      -DModule_ITKIOCSV:BOOL=ON
      -DModule_ITKIODCMTK:BOOL=ON
      -DModule_ITKIOGDCM:BOOL=ON
      -DModule_ITKIOGE:BOOL=ON
      -DModule_ITKIOGIPL:BOOL=ON
      -DModule_ITKIOHDF5:BOOL=ON # For DTIProcess
      -DModule_ITKIOIPL:BOOL=ON
      -DModule_ITKIOImageBase:BOOL=ON
      -DModule_ITKIOJPEG:BOOL=ON
      -DModule_ITKIOLSM:BOOL=ON
      -DModule_ITKIOMINC:BOOL=ON
      -DModule_ITKIOMRC:BOOL=ON
      -DModule_ITKIOMesh:BOOL=ON
      -DModule_ITKIOMeta:BOOL=ON
      -DModule_ITKIONIFTI:BOOL=ON
      -DModule_ITKIONRRD:BOOL=ON
      -DModule_ITKIOPNG:BOOL=ON
      -DModule_ITKIORAW:BOOL=ON
      -DModule_ITKIOSiemens:BOOL=ON
      -DModule_ITKIOSpatialObjects:BOOL=ON # For DTIProcess
      -DModule_ITKIOStimulate:BOOL=ON
      -DModule_ITKIOTIFF:BOOL=ON
      -DModule_ITKIOTransformBase:BOOL=ON
      -DModule_ITKIOTransformHDF5:BOOL=ON
      -DModule_ITKIOTransformInsightLegacy:BOOL=ON # For DTIProcess
      -DModule_ITKIOTransformMatlab:BOOL=ON
      -DModule_ITKIOVTK:BOOL=ON
      -DModule_ITKIOXML:BOOL=ON # For SlicerExecutionModel
      -DModule_ITKImageCompare:BOOL=ON
      -DModule_ITKMathematicalMorphology:BOOL=ON
      -DModule_ITKMetricsv4:BOOL=ON
      -DModule_ITKOptimizersv4:BOOL=ON
      -DModule_ITKRegionGrowing:BOOL=ON
      -DModule_ITKRegistrationCommon:BOOL=ON
      -DModule_ITKRegistrationMethodsv4:BOOL=ON
      -DModule_ITKTestKernel:BOOL=ON
      -DModule_ITKTransform:BOOL=ON
      -DModule_ITKVTK:BOOL=ON
      -DModule_MGHIO:BOOL=ON
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON
      -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      -DITK_WRAPPING:BOOL=OFF #${BUILD_SHARED_LIBS} ## HACK:  QUICK CHANGE
      # ZLIB
      -DITK_USE_SYSTEM_ZLIB:BOOL=ON
      -DZLIB_ROOT:PATH=${ZLIB_ROOT}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      -DITK_USE_FFTWD:BOOL=OFF
      -DITK_USE_FFTWF:BOOL=OFF
      ${ITK_VTK_OPTIONS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

# NOT NEEDED  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(ITK_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  if(NOT DEFINED ITK_VALGRIND_SUPPRESSIONS_FILE)
    set(ITK_VALGRIND_SUPPRESSIONS_FILE ${EP_SOURCE_DIR}/CMake/InsightValgrind.supp)
  endif()
  mark_as_superbuild(ITK_VALGRIND_SUPPRESSIONS_FILE:FILEPATH)

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  set(_lib_subdir lib)
  if(WIN32)
    set(_lib_subdir bin)
  endif()

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${ITK_DIR}/${_lib_subdir}/<CMAKE_CFG_INTDIR>)
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
    )

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS ITK_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
