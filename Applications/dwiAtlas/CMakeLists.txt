project(dwiAtlas)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

link_directories(${BOOST_LIBRARY_DIR} ${DTIIO_LIBRARY_DIR})

include_directories(${BOOST_INCLUDE_DIR} ${DWIProcess_SOURCE_DIR}/Library ${DWIProcess_SOURCE_DIR}/PrivateLibrary)

# Win32 has automatic linking of boost libraries so only add boost
# libaries on unix
# if(NOT WIN32 OR CYGWIN)
#   set(BOOST_LIBRARIES boost_program_options)
# endif(NOT WIN32 OR CYGWIN)

#find_package( ITK )
#if( ITK_FOUND )
#  include(${ITK_USE_FILE})
#else( ITK_FOUND )
#  message( FATAL_ERROR "Cannot build without ITK" )
#endif( ITK_FOUND )

#find_package(VTK)
#if(VTK_FOUND)
#  include(${VTK_USE_FILE})
#else(VTK_FOUND)
#  message(FATAL_ERROR
#  "VTK not found.  Please set VTK_DIR.")
#endif(VTK_FOUND)

#find_package(GenerateCLP REQUIRED)
#if(GenerateCLP_FOUND)
#  include(${GenerateCLP_USE_FILE})
#endif(GenerateCLP_FOUND)

set(dwiAtlas_SOURCE dwiAtlas.cxx)
GENERATECLP(dwiAtlas_SOURCE dwiAtlas.xml)

if(NOT EXECUTABLE_OUTPUT_PATH)
      set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
endif(NOT EXECUTABLE_OUTPUT_PATH)

add_executable(dwiAtlas dwiAtlas.cxx)

target_link_libraries(dwiAtlas ${ITK_LIBRARIES} ${ITK_LIBRARIES} ${ITK_LIBRARIES} ${ITK_LIBRARIES} DTIIO ${BOOST_LIBRARIES})


