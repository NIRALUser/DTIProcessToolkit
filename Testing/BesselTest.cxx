#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "wrapper.h"


//double args[] = {};
//double results[] = {};

int BesselTest(int argc, char* argv[])
{
  std::cout << "running Bessel Test" << std::endl;

  int ierr1, ierr2;
  for(double x = 0.0; x < 5; x += .1)
    {
    double unscaled =  besseli0(x,0,&ierr1);
    double scaled = besseli0(x,1,&ierr2);
    if(ierr1 || ierr2)
      {
      std::cout << "Internal error in bessel routine [FAILED]" << std::endl;
      return EXIT_FAILURE;
      }

    if(fabs(unscaled  - scaled * exp(x)) >= 1e-6)
      {
      std::cout << "Scaling of bessel function [FAILED]" << std::endl;
      return EXIT_FAILURE;
      }

    }
    
    
    std::cout << "Test [PASSED]. " << std::endl;
    return EXIT_SUCCESS;
}
