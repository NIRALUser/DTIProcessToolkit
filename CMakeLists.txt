PROJECT(DTIProcess)

NeuroLib_USE_ITK()
NeuroLib_USE_VTK()
NeuroLiB_USE_BOOST()

INCLUDE_DIRECTORIES(
${DTIProcess_SOURCE_DIR}/Library
${DTIProcess_SOURCE_DIR}/PrivateLibrary
)

# Workaround problems with cmake and fortran
ADD_SUBDIRECTORY(bessel)

ADD_SUBDIRECTORY(Library)
ADD_SUBDIRECTORY(PrivateLibrary)
ADD_SUBDIRECTORY(Applications)

ADD_SUBDIRECTORY(Testing)
