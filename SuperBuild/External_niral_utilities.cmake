
set(proj niral_utilities)
# Set dependency list
set(${proj}_DEPENDENCIES "zlib")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)


# Sanity checks
if(DEFINED niral_utilities_DIR AND NOT EXISTS ${niral_utilities_DIR})
  message(FATAL_ERROR "niral_utilities_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED niral_utilities_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/NIRALUser/niral_utilities.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "release" # 
    QUIET
    )

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  ExternalProject_Add(${proj}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    ${${proj}_EP_ARGS}
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
##  niral_utilities Flags
    -DCOMPILE_CONVERTITKFORMATS:BOOL=OFF
    -DCOMPILE_CROPTOOLS:BOOL=${BUILD_CropDTI}
    -DCOMPILE_CURVECOMPARE:BOOL=OFF
    -DCOMPILE_DWI_NIFTINRRDCONVERSION:BOOL=OFF
    -DCOMPILE_IMAGEMATH:BOOL=${COMPILE_IMAGEMATH}
    -DCOMPILE_IMAGESTAT:BOOL=OFF
    -DCOMPILE_MULTIATLASSEG:BOOL=OFF
    -DCOMPILE_POLYDATAMERGE:BOOL=${BUILD_PolyDataMerge}
    -DCOMPILE_POLYDATATRANSFORM:BOOL=${BUILD_PolyDataTransform}
    -DCOMPILE_Unsupported:BOOL=OFF
    -DUSE_SYSTEM_SlicerExecutionModel:BOOL=ON
    -DUSE_SYSTEM_ITK:BOOL=ON
    -DUSE_SYSTEM_VTK:BOOL=ON
    -DITK_VERSION_MAJOR:STRING=${ITK_VERSION_MAJOR}
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
## INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )

# NOT NEEDED  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(${extProjName}_DIR ${EXTERNAL_BINARY_DIRECTORY}/${proj}-install/lib/CMake/${proj})

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS niral_utilities_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
