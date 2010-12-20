ADD_EXECUTABLE(nrrdvtk nrrdvtk.cxx)
TARGET_LINK_LIBRARIES(nrrdvtk ITKIO)

ADD_EXECUTABLE(dwiprocess dwiprocess.cxx)
TARGET_LINK_LIBRARIES(dwiprocess ITKIO)

ADD_EXECUTABLE(labeltable labeltable.cxx)
TARGET_LINK_LIBRARIES(labeltable ITKIO)

ADD_EXECUTABLE(hfieldtodeformation  hfieldtodeformation.cxx)
TARGET_LINK_LIBRARIES(hfieldtodeformation TensorOperations DTIIO ITKIO)

ADD_EXECUTABLE(tio tio.cxx)
TARGET_LINK_LIBRARIES(tio ITKIO)

ADD_EXECUTABLE(morph morph.cxx)
TARGET_LINK_LIBRARIES(morph ITKIO)

ADD_SUBDIRECTORY(dwiAtlas)
ADD_DEPENDENCIES(dwiAtlas dwiprocess DTIIO)	


set(dtiprocess_SRC dtiprocess.cxx)
GENERATECLP(dtiprocess_SRC dtiprocess.xml)
ADD_EXECUTABLE(dtiprocess ${dtiprocess_SRC})
TARGET_LINK_LIBRARIES(dtiprocess TensorOperations DTIIO ${ITK_LIBRARIES})

set(dtiestim_SRC dtiestim.cxx)
GENERATECLP(dtiestim_SRC dtiestim.xml)
ADD_EXECUTABLE(dtiestim ${dtiestim_SRC})
TARGET_LINK_LIBRARIES(dtiestim DWIPrivateLibrary ${ITK_LIBRARIES} cephes)

set(scalartransform_SRC scalartransform.cxx)
GENERATECLP(scalartransform_SRC scalartransform.xml)
ADD_EXECUTABLE(scalartransform ${scalartransform_SRC})
TARGET_LINK_LIBRARIES(scalartransform DTIIO ITKIO )

set(fiberprocess_SRC fiberprocess.cxx)
GENERATECLP(fiberprocess_SRC fiberprocess.xml)
ADD_EXECUTABLE(fiberprocess ${fiberprocess_SRC})
TARGET_LINK_LIBRARIES(fiberprocess DTIIO ITKIO ITKNumerics )

set(dtiaverage_SRC dtiaverage.cxx)
GENERATECLP(dtiaverage_SRC dtiaverage.xml)
ADD_EXECUTABLE(dtiaverage ${dtiaverage_SRC})
TARGET_LINK_LIBRARIES(dtiaverage TensorOperations DTIIO )

set(fibertrack_SRC fibertrack.cxx)
GENERATECLP(fibertrack fibertrack.xml)
ADD_EXECUTABLE(fibertrack ${fibertrack_SRC})
TARGET_LINK_LIBRARIES(fibertrack DTIIO TensorOperations ITKIO ITKSpatialObject )

set(fiberstats_SRC fiberstats.cxx)
GENERATECLP(fiberstats_SRC fiberstats.xml)
ADD_EXECUTABLE(fiberstats fiberstats.cxx)
TARGET_LINK_LIBRARIES(fiberstats DTIIO )

set(maxcurvature_SRC maxcurvature.cxx)
GENERATECLP(maxcurvature_SRC maxcurvature.xml)
ADD_EXECUTABLE(maxcurvature ${maxcurvature_SRC})
TARGET_LINK_LIBRARIES(maxcurvature ITKIO ITKNumerics )


INSTALL(TARGETS dtiprocess 
                dtiestim
                scalartransform
                #                fiberprocess
                dtiaverage
                fibertrack
                fiberstats
                maxcurvature
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib/static)

