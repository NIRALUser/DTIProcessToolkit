
ADD_LIBRARY(DWIPrivateLibrary wrapper.cxx)

TARGET_LINK_LIBRARIES(DWIPrivateLibrary ${BESSEL_LIB} )

ADD_LIBRARY(TensorOperations tensorscalars.cxx tensordeformation.cxx)
ADD_LIBRARY(DTIIO tensorio.cxx fiberio.cxx deformationfieldio.cxx)
TARGET_LINK_LIBRARIES(DTIIO vtkIO ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(TensorOperations vtkIO ${ITK_LIBRARIES})

set( libraries_targets DWIPrivateLibrary TensorOperations DTIIO )

set_target_properties(${libraries_targets} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${ARCHIVE_OUTPUT_DIRECTORY}"
    )

install(TARGETS ${libraries_targets}
    LIBRARY DESTINATION "${Slicer_INSTALL_ROOT}${Slicer_BUNDLE_EXTENSIONS_LOCATION}${Slicer_LIB_DIR}" COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION "${Slicer_INSTALL_ROOT}${Slicer_BUNDLE_EXTENSIONS_LOCATION}${Slicer_LIB_DIR}" COMPONENT RuntimeLibraries
    )
