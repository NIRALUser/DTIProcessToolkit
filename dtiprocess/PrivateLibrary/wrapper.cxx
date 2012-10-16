#include "wrapper.h"

#if 0

extern "C" void zbesi_(double *zr,
                   double *zi,
                   double *fnu,
                   int *kode,
                   int *n,
                   double *cyr,
                   double *cyi,
                   int *nz,
                   int *nerr);

// If scaling is non-zero then the result * exp(-abs(x))
// gives the bessel function
// This function returns exp(-abs(x))*I(0,x);
double besseli0(double x, int scaling, int *ierr_)
{
  double zero = 0.0;
  int n = 1;
  int kode = !scaling ? 1 : 2;
  double cyr, cyi;
  int nz, ierr;

  zbesi_(&x,&zero,&zero,&kode,&n,&cyr,&cyi,&nz,&ierr);
  *ierr_ = ierr;

  return cyr;
}

#if 0
double besseli1(double x, int scaling, int *ierr_)
{
  double zero = 0.0;
  double one = 1.0;
  int n = 1;
  int kode = !scaling ? 1 : 2;
  double cyr, cyi;
  int nz, ierr;

  zbesi_(&x,&zero,&one,&kode,&n,&cyr,&cyi,&nz,&ierr);
  *ierr_ = ierr;

  return cyr;
}
#endif
#endif
