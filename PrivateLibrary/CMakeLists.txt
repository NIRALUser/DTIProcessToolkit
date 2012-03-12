
ADD_LIBRARY(DWIPrivateLibrary wrapper.cxx)

TARGET_LINK_LIBRARIES(DWIPrivateLibrary ${BESSEL_LIB} )

ADD_LIBRARY(TensorOperations tensorscalars.cxx tensordeformation.cxx)
ADD_LIBRARY(DTIIO tensorio.cxx fiberio.cxx deformationfieldio.cxx)
TARGET_LINK_LIBRARIES(DTIIO vtkIO ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(TensorOperations vtkIO ${ITK_LIBRARIES})

INSTALL(TARGETS DWIPrivateLibrary TensorOperations DTIIO
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib)
