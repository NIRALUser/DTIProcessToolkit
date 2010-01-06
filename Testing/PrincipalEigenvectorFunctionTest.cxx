/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/14 15:54:58 $
  Version:   $Revision: 1.1 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <cstdlib>

#include <itkDiffusionTensor3D.h>


#include "itkTensorPrincipalEigenvectorImageFilter.h"


int PrincipalEigenvectorFunctionTest(int argc, char* argv[])
{
  std::cout << "running Principal Eigenvector test" << std::endl;
  
  typedef itk::DiffusionTensor3D<double> TensorType;
  TensorType a(0.0);
  a[0] = 3; a[3] = 1.5; a[5] = 1.0;

  typedef itk::Functor::TensorPrincipalEigenvectorFunction<TensorType, double>  PrincipalEigenvectorFunctorType;

  PrincipalEigenvectorFunctorType findprineig;

  typedef PrincipalEigenvectorFunctorType::PixelType OutputType;
  OutputType x = findprineig(a);

  OutputType correctx;
  correctx[0] = 1; correctx[1] = 0; correctx[2] = 0;

  for(unsigned int i = 0; i < 3; i++)
  {
    if(fabs(x[i] - correctx[i]) > 1e-6)
    {
      std::cout << "Principal eigenvector incorrect [FAILED]" << std::endl;
      std::cout << "Correct: " << correctx << std::endl;
      std::cout << "Found  : " << x << std::endl;
      return EXIT_FAILURE;
    }
  }


  std::cout << "Test [PASSED]." << std::endl;
  return EXIT_SUCCESS;
}
