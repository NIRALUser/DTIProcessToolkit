

ADD_LIBRARY(DWIPrivateLibrary ${STATIC_LIB} wrapper.cxx)

TARGET_LINK_LIBRARIES(DWIPrivateLibrary ${BESSEL_LIB} )

ADD_LIBRARY(TensorOperations ${STATIC_LIB} tensorscalars.cxx tensordeformation.cxx)
ADD_LIBRARY(DTIIO ${STATIC_LIB} tensorio.cxx fiberio.cxx deformationfieldio.cxx)
TARGET_LINK_LIBRARIES(DTIIO vtkIO ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(TensorOperations vtkIO ${ITK_LIBRARIES})

set( libraries_targets DWIPrivateLibrary TensorOperations DTIIO )

set_target_properties(${libraries_targets} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    )

install(TARGETS ${libraries_targets}
    RUNTIME DESTINATION "${INSTALL_RUNTIME_DESTINATION}" COMPONENT RuntimeLibraries
    LIBRARY DESTINATION "${INSTALL_LIBRARY_DESTINATION}" COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION "${INSTALL_ARCHIVE_DESTINATION}" COMPONENT RuntimeLibraries
    )

