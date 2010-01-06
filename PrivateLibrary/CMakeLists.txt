
ADD_LIBRARY(DWIPrivateLibrary wrapper.cxx)

TARGET_LINK_LIBRARIES(DWIPrivateLibrary ${BESSEL_LIB} )

ADD_LIBRARY(TensorOperations SHARED tensorscalars.cxx tensordeformation.cxx)
ADD_LIBRARY(DTIIO SHARED tensorio.cxx fiberio.cxx deformationfieldio.cxx)
TARGET_LINK_LIBRARIES(DTIIO vtkIO ITKIO)
TARGET_LINK_LIBRARIES(TensorOperations vtkIO ITKIO)

INSTALL(TARGETS DWIPrivateLibrary TensorOperations DTIIO 
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib)
