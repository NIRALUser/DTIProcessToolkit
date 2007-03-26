NeuroLib_USE_ITK()
NeuroLib_USE_BOOST()

ADD_EXECUTABLE(dtiestim dtiestim.cxx) # wrapper.cxx)
ADD_EXECUTABLE(dtiprocess dtiprocess.cxx)
ADD_EXECUTABLE(nrrdvtk nrrdvtk.cxx)
ADD_EXECUTABLE(dtiaverage dtiaverage.cxx)
ADD_EXECUTABLE(maxcurvature maxcurvature.cxx)
ADD_EXECUTABLE(fiberprocess fiberprocess.cxx)

TARGET_LINK_LIBRARIES(dtiestim ITKCommon ITKNumerics ITKNrrdIO ITKIO
ITKStatistics itkvnl ${BOOST_LIBRARIES} bessel gfortran)
TARGET_LINK_LIBRARIES(fiberprocess ITKCommon ITKNumerics ITKNrrdIO ITKIO
ITKStatistics ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(dtiprocess ITKCommon ITKNumerics ITKNrrdIO ITKIO
ITKStatistics ${BOOST_LIBRARIES})
TARGET_LINK_LIBRARIES(nrrdvtk ITKCommon ITKNumerics ITKNrrdIO ITKIO
ITKStatistics)
TARGET_LINK_LIBRARIES(dtiaverage ITKCommon ITKNumerics ITKNrrdIO ITKIO
ITKStatistics)
TARGET_LINK_LIBRARIES(maxcurvature ITKCommon ITKNumerics ITKNrrdIO ITKIO
ITKStatistics)

