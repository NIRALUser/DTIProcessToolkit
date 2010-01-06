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

#include "transforms.h"

#include <itkAffineTransform.h>

int ReadITKAffineTest(int argc, char* argv[])
{
  const std::string floatfile(argv[1]);
  const std::string doublefile(argv[2]);

  itk::AffineTransform<float, 3>::Pointer floatfromfloat =
    itk::AffineTransform<float, 3>::New();
  itk::AffineTransform<double, 3>::Pointer doublefromdouble =
    itk::AffineTransform<double, 3>::New();

  itk::AffineTransform<float, 3>::Pointer floatfromdouble =
    itk::AffineTransform<float, 3>::New();
  itk::AffineTransform<double, 3>::Pointer doublefromfloat =
    itk::AffineTransform<double, 3>::New();

  std::cout << "Reading float transform into float transform" << std::endl;
  floatfromfloat = readITKAffine<float, 3>(floatfile);
  std::cout << floatfromfloat->GetParameters() << std::endl;
  std::cout << floatfromfloat->GetFixedParameters() << std::endl;

  std::cout << "Reading double transform into double transform" << std::endl;
  doublefromdouble = readITKAffine<double, 3>(doublefile);
  std::cout << doublefromdouble->GetParameters() << std::endl;
  std::cout << doublefromdouble->GetFixedParameters() << std::endl;

  std::cout << "Reading double transform into float transform" << std::endl;
  floatfromdouble = readITKAffine<float, 3>(doublefile);
  std::cout << floatfromdouble->GetParameters() << std::endl;
  std::cout << floatfromdouble->GetFixedParameters() << std::endl;

  std::cout << "Reading float transform into double transform" << std::endl;
  doublefromfloat = readITKAffine<double, 3>(floatfile);
  std::cout << doublefromfloat->GetParameters() << std::endl;
  std::cout << doublefromfloat->GetFixedParameters() << std::endl;

  return 0;
}
