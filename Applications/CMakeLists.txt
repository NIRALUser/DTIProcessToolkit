FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp)

INCLUDE_DIRECTORIES(
${BOOST_INCLUDE_DIR}
)

FIND_PATH(BOOST_LIBRARY_DIR boost_date_time-vc71-mt.lib)
LINK_DIRECTORIES(${BOOST_LIBRARY_DIR})

ADD_EXECUTABLE(dtiestim dtiestim.cxx)
ADD_EXECUTABLE(dtiprocess dtiprocess.cxx)
ADD_EXECUTABLE(nrrdvtk nrrdvtk.cxx)
ADD_EXECUTABLE(dtiaverage dtiaverage.cxx)
ADD_EXECUTABLE(maxcurvature maxcurvature.cxx)
ADD_EXECUTABLE(dtipgaaverage dtipgaaverage.cxx)
ADD_EXECUTABLE(fiberprocess fiberprocess.cxx)
ADD_EXECUTABLE(roiprocess roiprocess.cxx)
ADD_EXECUTABLE(tensorfeature tensorfeature.cxx)
ADD_EXECUTABLE(dwiprocess dwiprocess.cxx)
ADD_EXECUTABLE(labeltable labeltable.cxx)
ADD_EXECUTABLE(hfieldtodeformation  hfieldtodeformation.cxx)
ADD_EXECUTABLE(tio tio.cxx)
ADD_EXECUTABLE(fibertrack fibertrack.cxx)
ADD_EXECUTABLE(morph morph.cxx)
ADD_EXECUTABLE(fiberstats fiberstats.cxx)
#ADD_EXECUTABLE(scalartransform scalartransform.cxx)

# Win32 has automatic linking of boost libraries so only add boost
# libaries on unix
IF(NOT WIN32 OR CYGWIN)
  SET(BOOST_LIBRARIES boost_program_options)
ENDIF(NOT WIN32 OR CYGWIN)

# SET(DWIPROCESS_DEFAULT_LIBRARIES ITKStatistics ITKNumerics ITKIO ITKCommon itkvnl DWIPrivateLibrary)

TARGET_LINK_LIBRARIES(dtiestim DWIPrivateLibrary ITKIO ITKNumerics ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(fiberprocess DTIIO ITKIO ITKNumerics ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(roiprocess ITKIO ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(dtiprocess TensorOperations DTIIO ITKIO ITKNumerics ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(nrrdvtk ITKIO)
TARGET_LINK_LIBRARIES(dtiaverage TensorOperations DTIIO ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(maxcurvature ITKIO ITKNumerics ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(dtipgaaverage TensorOperations DTIIO)
TARGET_LINK_LIBRARIES(tensorfeature ITKIO)
TARGET_LINK_LIBRARIES(dwiprocess ITKIO)
TARGET_LINK_LIBRARIES(labeltable ITKIO)
TARGET_LINK_LIBRARIES(hfieldtodeformation TensorOperations DTIIO ITKIO)
TARGET_LINK_LIBRARIES(tio ITKIO)
TARGET_LINK_LIBRARIES(fibertrack DTIIO TensorOperations ITKIO
ITKSpatialObject ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(morph ITKIO)
TARGET_LINK_LIBRARIES(fiberstats DTIIO ${BOOST_LIBRARIES})
#TARGET_LINK_LIBRARIES(scalartransform DTIIO ITKIO ${BOOST_LIBRARIES})
