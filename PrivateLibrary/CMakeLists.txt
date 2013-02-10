
add_library(TensorOperations tensorscalars.cxx tensordeformation.cxx)
add_library(DTIIO tensorio.cxx fiberio.cxx deformationfieldio.cxx)
target_link_libraries(DTIIO vtkIO ${ITK_LIBRARIES})
target_link_libraries(TensorOperations vtkIO ${ITK_LIBRARIES})

install(TARGETS TensorOperations DTIIO
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib)
