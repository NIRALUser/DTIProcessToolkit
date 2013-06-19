
set(cephes_SRC i0.c)
add_library(cephes ${cephes_SRC})


set( libraries_targets cephes )


set_target_properties(${libraries_targets} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    )

install(TARGETS ${libraries_targets}
    LIBRARY DESTINATION "${Slicer_INSTALL_ROOT}${Slicer_BUNDLE_EXTENSIONS_LOCATION}${Slicer_LIB_DIR}" COMPONENT RuntimeLibraries
    )

