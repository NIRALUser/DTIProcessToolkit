/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorApparentDiffusionCoefficientImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorApparentDiffusionCoefficientImageFilter_h
#define __itkTensorApparentDiffusionCoefficientImageFilter_h

#include <itkDiffusionTensor3D.h>

namespace itk
{

namespace Functor 
{

/**
 *  Compute the apparent diffusion coefficient of a tensor along a specific direction.
 */
template <typename TInput, typename TGradientDirection, typename TOutput>
class TensorApparentDiffusionCoefficient
{
public:
  TensorApparentDiffusionCoefficient() {}
  ~TensorApparentDiffusionCoefficient() {}

  bool operator!=( const TensorApparentDiffusionCoefficient & ) const
  {
    return false;
  }
  bool operator==( const TensorApparentDiffusionCoefficient & other ) const
  {
    return !(*this != other);
  }

  TOutput operator()(const TInput & d, const TGradientDirection & direction)
  {
    double result = 0.0;
    for(unsigned int i = 0; i < 3; ++i)
      for(unsigned int j = 0; j < 3; ++j)
        result += direction[i]*d(i,j)*direction[j];
    return result;
  }
  
};

} // end namespace Functor

} // end namespace itk
#endif
