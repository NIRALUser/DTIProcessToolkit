#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "cephes/cephes.h"


// bessel function replaced with cephes #include "wrapper.h"
//double args[] = {};
//double results[] = {};

int BesselTest(int argc, char* argv[])
{
  std::cout << "running Bessel Test" << std::endl;

  for(double x = 0.0; x < 5; x += .1)
    {
    //int ierr1, ierr2;
    // bessel function replaced with cephes double unscaled =  besseli0(x,0,&ierr1);
    // bessel function replaced with cephes double scaled = besseli0(x,1,&ierr2);
    // if(ierr1 || ierr2)
    //  {
    //  std::cout << "Internal error in bessel routine [FAILED]" << std::endl;
    //  return EXIT_FAILURE;
    //  }

    const double unscaled =  i0(x);
    const double scaled = i0e(x);
    if(fabs(unscaled  - scaled * exp(x)) >= 1e-6)
      {
      std::cout << "Scaling of bessel function [FAILED]" << std::endl;
      return EXIT_FAILURE;
      }

    }
    std::cout << "Test [PASSED]. " << std::endl;
    return EXIT_SUCCESS;
}
