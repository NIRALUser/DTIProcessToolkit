#include "wrapper.h"

extern "C" void zbesi_(double *zr,
                   double *zi,
                   double *fnu,
                   int *kode,
                   int *n,
                   double *cyr,
                   double *cyi,
                   int *nz,
                   int *nerr);

// If scaling is non-zero then the result * exp(x)
// gives the bessel function
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

