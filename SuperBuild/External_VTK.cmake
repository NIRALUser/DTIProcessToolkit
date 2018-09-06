if( NOT EXTERNAL_SOURCE_DIRECTORY )
  set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
endif()

# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

## External_${extProjName}.cmake files can be recurisvely included,
## and cmake variables are global, so when including sub projects it
## is important make the extProjName and proj variables
## appear to stay constant in one of these files.
## Store global variables before overwriting (then restore at end of this file.)
ProjectDependancyPush(CACHED_extProjName ${extProjName})
ProjectDependancyPush(CACHED_proj ${proj})

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName VTK) #The find_package known name
set(proj        VTK) #This local name


#Setting VTK_VERSION_MAJOR to its default value if it has not been set before
# set(VTK_VERSION_MAJOR 5 CACHE STRING "Choose the expected VTK major version to build Slicer (5 or 6).")
# # Set the possible values of VTK major version for cmake-gui
# set_property(CACHE VTK_VERSION_MAJOR PROPERTY STRINGS "5" "6")
# if(NOT "${VTK_VERSION_MAJOR}" STREQUAL "5" AND NOT "${VTK_VERSION_MAJOR}" STREQUAL "6")
#   set(VTK_VERSION_MAJOR 5 CACHE STRING "Choose the expected VTK major version to build Slicer (5 or 6)." FORCE)
#   message(WARNING "Setting VTK_VERSION_MAJOR to '5' as an valid value was specified.")
# endif()

# set(USE_VTKv5 ON)
# set(USE_VTKv6 OFF)
# if(${VTK_VERSION_MAJOR} STREQUAL "6")
#   set(USE_VTKv5 OFF)
#   set(USE_VTKv6 ON)
# endif()

# if(USE_VTKv6)
#   set(${extProjName}_REQUIRED_VERSION "6.1")  #If a required version is necessary, then set this, else leave blank
# else()
#   set(${extProjName}_REQUIRED_VERSION "5.10")  #If a required version is necessary, then set this, else leave blank
# endif()

# Sanity checks
#if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
#  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
#endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")
set(${PROJECT_NAME}_USE_PYTHONQT OFF)
if (${PROJECT_NAME}_USE_PYTHONQT)
  list(APPEND ${proj}_DEPENDENCIES python)
endif()

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
      -DVTK_REQUIRED_OBJCXX_FLAGS:STRING="")
  endif()

  ### --- Project specific additions here
  set(VTK_WRAP_TCL OFF)
  set(VTK_WRAP_PYTHON OFF)

  if (${PROJECT_NAME}_USE_PYTHONQT)
    set(VTK_WRAP_PYTHON ON)
  endif()

  set(VTK_PYTHON_ARGS
      -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
      -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
      -DPYTHON_LIBRARIES:FILEPATH=${PYTHON_LIBRARIES}
      )
  if(${PROJECT_NAME}_USE_PYTHONQT)
    list(APPEND VTK_PYTHON_ARGS
      -DVTK_INSTALL_PYTHON_USING_CMAKE:BOOL=ON
      )
  endif()

  set(VTK_QT_ARGS)
  if(${PRIMARY_PROJECT_NAME}_USE_QT)
    # if(USE_VTKv6)
      set(VTK_QT_ARGS
        -DModule_vtkGUISupportQt:BOOL=ON
        )
    # else()
    #   set(VTK_QT_ARGS)
    # endif()
    if(NOT APPLE)
      list(APPEND VTK_QT_ARGS
        #-DDESIRED_QT_VERSION:STRING=4 # Unused
        -DVTK_USE_GUISUPPORT:BOOL=ON
        -DVTK_USE_QVTK_QTOPENGL:BOOL=ON
        -DVTK_USE_QT:BOOL=ON
        -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        )
    else()
      list(APPEND VTK_QT_ARGS
        -DVTK_USE_CARBON:BOOL=OFF
        # Default to Cocoa, VTK/CMakeLists.txt will enable Carbon and disable cocoa if needed
        -DVTK_USE_COCOA:BOOL=ON
        -DVTK_USE_X:BOOL=OFF
        #-DVTK_USE_RPATH:BOOL=ON # Unused
        #-DDESIRED_QT_VERSION:STRING=4 # Unused
        -DVTK_USE_GUISUPPORT:BOOL=ON
        -DVTK_USE_QVTK_QTOPENGL:BOOL=ON
        -DVTK_USE_QT:BOOL=ON
        -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        )
    endif()
    find_package(Qt4 REQUIRED)
  else()
    set(VTK_QT_ARGS
        -DVTK_USE_GUISUPPORT:BOOL=OFF
        -DVTK_USE_QT:BOOL=OFF
        )
  endif()

  # Disable Tk when Python wrapping is enabled
  if (${PROJECT_NAME}_USE_PYTHONQT)
    list(APPEND VTK_QT_ARGS -DVTK_USE_TK:BOOL=OFF)
  endif()

  set(slicer_TCL_LIB)
  set(slicer_TK_LIB)
  set(slicer_TCLSH)
  set(VTK_TCL_ARGS)
  if(VTK_WRAP_TCL)
    if(WIN32)
      set(slicer_TCL_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/tcl84.lib)
      set(slicer_TK_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/tk84.lib)
      set(slicer_TCLSH ${CMAKE_BINARY_DIR}/tcl-build/bin/tclsh.exe)
    elseif(APPLE)
      set(slicer_TCL_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtcl8.4.dylib)
      set(slicer_TK_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtk8.4.dylib)
      set(slicer_TCLSH ${CMAKE_BINARY_DIR}/tcl-build/bin/tclsh84)
    else()
      set(slicer_TCL_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtcl8.4.so)
      set(slicer_TK_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtk8.4.so)
      set(slicer_TCLSH ${CMAKE_BINARY_DIR}/tcl-build/bin/tclsh84)
    endif()
    set(VTK_TCL_ARGS
      -DTCL_INCLUDE_PATH:PATH=${CMAKE_BINARY_DIR}/tcl-build/include
      -DTK_INCLUDE_PATH:PATH=${CMAKE_BINARY_DIR}/tcl-build/include
      -DTCL_LIBRARY:FILEPATH=${slicer_TCL_LIB}
      -DTK_LIBRARY:FILEPATH=${slicer_TK_LIB}
      -DTCL_TCLSH:FILEPATH=${slicer_TCLSH}
      )
  endif()

  set(VTK_BUILD_STEP "")
  if(UNIX)
    configure_file(SuperBuild/External_VTK_build_step.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/External_VTK_build_step.cmake
      @ONLY)
    set(VTK_BUILD_STEP ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/External_VTK_build_step.cmake)
  endif()

  set(${proj}_CMAKE_OPTIONS
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${${PROJECT_NAME}_USE_VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=OFF
      -DVTK_WRAP_TCL:BOOL=${VTK_WRAP_TCL}
      #-DVTK_USE_RPATH:BOOL=ON # Unused
      ${VTK_TCL_ARGS}
      -DVTK_WRAP_PYTHON:BOOL=${VTK_WRAP_PYTHON}
      -DVTK_INSTALL_LIB_DIR:PATH=${${PROJECT_NAME}_INSTALL_LIB_DIR}
      ${VTK_PYTHON_ARGS}
      ${VTK_QT_ARGS}
      ${VTK_MAC_ARGS}
    )
  ### --- End Project specific additions
  # if(USE_VTKv6)
  set(${proj}_GIT_TAG "v7.0.0")
  set(${proj}_REPOSITORY ${git_protocol}://vtk.org/VTK.git)
  # else()
  #   set(${proj}_REPOSITORY ${git_protocol}://github.com/BRAINSia/VTK.git)
  #   set(${proj}_GIT_TAG "FixClangFailure_VTK5.10_release")
  # endif()
  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${EXTERNAL_SOURCE_DIRECTORY}/${proj}
    BINARY_DIR ${proj}-build
    BUILD_COMMAND ${VTK_BUILD_STEP}
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      ${${proj}_CMAKE_OPTIONS}
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )


  # set(VTKPatchScript ${CMAKE_CURRENT_LIST_DIR}/External_VTK_patch.cmake)
  # ExternalProject_Add_Step(${proj} VTKPatch
  #   COMMENT "get rid of obsolete C/CXX flags"
  #   DEPENDEES download
  #   DEPENDERS configure
  #   COMMAND ${CMAKE_COMMAND}
  #   -DVTKSource=<SOURCE_DIR>
  #   -DUSE_VTKv6=${USE_VTKv6}
  #   -P ${VTKPatchScript}
  #   )

# if(USE_VTKv6)
#   set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/cmake/vtk-6.1)
# else()
#   set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/vtk-5.10)
# endif()

set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/cmake/vtk-7.0)

else()
  if(${USE_SYSTEM_${extProjName}})
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
    find_package(${extProjName} REQUIRED)
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
