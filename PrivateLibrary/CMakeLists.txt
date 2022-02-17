

add_library(TensorOperations ${STATIC_LIB} tensorscalars.cxx tensordeformation.cxx)
add_library(DTIIO ${STATIC_LIB} tensorio.cxx fiberio.cxx deformationfieldio.cxx)
target_link_libraries(DTIIO ${VTK_LIBRARIES} ${ITK_LIBRARIES})
target_link_libraries(TensorOperations ${VTK_LIBRARIES} ${ITK_LIBRARIES})

set( libraries_targets TensorOperations DTIIO )

set_target_properties(${libraries_targets} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    )

install(TARGETS ${libraries_targets}
    RUNTIME DESTINATION "${INSTALL_RUNTIME_DESTINATION}" COMPONENT RuntimeLibraries
    LIBRARY DESTINATION "${INSTALL_LIBRARY_DESTINATION}" COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION "${INSTALL_ARCHIVE_DESTINATION}" COMPONENT RuntimeLibraries
    )

