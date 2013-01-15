/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionWeightedImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionWeightedImageFilter_txx
#define __itkDiffusionTensor3DReconstructionWeightedImageFilter_txx

#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_svd.h"

namespace itk {

template< class TGradientImagePixelType,
          class TTensorPrecision >
DiffusionTensor3DReconstructionWeightedImageFilter< TGradientImagePixelType, 
                                                    TTensorPrecision >
::DiffusionTensor3DReconstructionWeightedImageFilter() :  m_NumberOfIterations(1)
{
}

template< class TGradientImagePixelType, class TTensorPrecision >
vnl_vector<TTensorPrecision>
DiffusionTensor3DReconstructionWeightedImageFilter< TGradientImagePixelType, 
                                                    TTensorPrecision >
::EstimateTensor(const vnl_vector<TTensorPrecision>& S) const
{
  // setup log signals
  vnl_vector< TTensorPrecision > B(this->m_NumberOfGradientDirections);
  for(unsigned int i = 0; i < S.size(); ++i)
  {
    if(S[i] == 0)
      B[i] = 0;
    else
      B[i] = log(S[i]);
  }

  vnl_vector<double> prevestimate = this->m_TensorBasis * B;
  for(unsigned int iter = 0; iter < this->m_NumberOfIterations; ++iter)
  {
    vnl_vector<double> phi(this->m_NumberOfGradientDirections);
    phi = this->m_BMatrix * prevestimate;
    
    for(unsigned int i = 0; i < this->m_NumberOfGradientDirections; ++i)
    {
      phi[i] = exp(phi[i]);
      phi[i] *= phi[i];
    }           
    
    vnl_diag_matrix<double> W2(phi);
    
    prevestimate = 
      vnl_svd<double>(this->m_BMatrix.transpose() * W2 * this->m_BMatrix).solve(this->m_BMatrix.transpose() * W2 * B);
    
  }

  return prevestimate;
}

}
#endif
