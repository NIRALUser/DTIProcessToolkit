PROJECT(bessel Fortran)

FILE(WRITE ${NeuroLib_BINARY_DIR}/THE.mod.proxy "")
FILE(WRITE ${NeuroLib_BINARY_DIR}/WITH.mod.proxy "")
FILE(WRITE ${NeuroLib_BINARY_DIR}/FORWARD.mod.proxy "")

ADD_LIBRARY(bessel d1mach.f  xzabs.f  zasyi.f  zbuni.f  zlog.f   zs1s2.f  zuchk.f	zunik.f
dgamln.f  zabs.f   zbesi.f  zdiv.f   zmlri.f  zseri.f  zunhj.f	zuoik.f
i1mach.f  zacai.f  zbinu.f  zexp.f   zmlt.f   zshch.f  zuni1.f	zwrsk.f
xerror.f  zairy.f  zbknu.f  zkscl.f  zrati.f  zsqrt.f  zuni2.f
)

