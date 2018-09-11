if( NOT EXTERNAL_SOURCE_DIRECTORY )
  set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
endif()
if( NOT EXTERNAL_BINARY_DIRECTORY )
  set( EXTERNAL_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
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
set(extProjName Boost) #The find_package known name
set(proj      Boost) #This local name

set(${proj}_DEPENDENCIES "")


if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  set(Boost_Install_Dir ${CMAKE_CURRENT_BINARY_DIR}/boost-install)
  set(Boost_Configure_Script ${CMAKE_CURRENT_LIST_DIR}/configureboost.cmake)
  set(Boost_Build_Script ${CMAKE_CURRENT_LIST_DIR}/buildboost.cmake)

  ExternalProject_add(Boost
    SVN_REPOSITORY http://svn.boost.org/svn/boost/trunk
    SVN_REVISION -r "77936"
  #  URL http://sourceforge.net/projects/boost/files/boost/1.49.0/boost_1_49_0.tar.gz
  #  URL_MD5 e0defc8c818e4f1c5bbb29d0292b76ca
    SOURCE_DIR Boost
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
    -DBUILD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/Boost
    -DBOOST_INSTALL_DIR:PATH=${Boost_Install_Dir}
    -P ${Boost_Configure_Script}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_COMMAND}
    -DBUILD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/Boost
    -DBOOST_INSTALL_DIR:PATH=${Boost_Install_Dir} -P ${Boost_Build_Script}
  )
  set(BOOST_ROOT ${CMAKE_CURRENT_BINARY_DIR}/boost-install)
  set(BOOST_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/boost-install/include)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} )
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS BOOST_ROOT:PATH BOOST_INCLUDE_DIR:PATH )
_expand_external_project_vars()
set(COMMON_EXTERNAL_PROJECT_ARGS ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
