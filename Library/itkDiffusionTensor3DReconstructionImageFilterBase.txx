/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.3 $

  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDiffusionTensor3DReconstructionImageFilterBase_txx
#define __itkDiffusionTensor3DReconstructionImageFilterBase_txx

#include "itkDiffusionTensor3DReconstructionImageFilterBase.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "vnl/vnl_vector.h"

namespace itk {

template< class TGradientImagePixelType, class TTensorPrecision >
DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                               TTensorPrecision >
::DiffusionTensor3DReconstructionImageFilterBase()
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 ); 
  m_NumberOfGradientDirections = 0;
  m_Threshold = NumericTraits< GradientPixelType >::min();
  m_GradientDirectionContainer = NULL;
  m_BValue = 1.0;
  m_EstimateBaseline = false;
}

template< class TGradientImagePixelType, class TTensorPrecision >
void DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                    TTensorPrecision >
::SetEstimateBaseline(bool eb)
{
  if( eb == m_EstimateBaseline )
    return;
  m_EstimateBaseline = eb;
  if(m_EstimateBaseline)
  {
    typename ScalarImageType::Pointer output = ScalarImageType::New();
    this->ProcessObject::SetNumberOfRequiredOutputs(2);
    this->ProcessObject::SetNthOutput(1, output.GetPointer());
  }
  else
  {
    this->ProcessObject::SetNthOutput(1, NULL);
    this->ProcessObject::SetNumberOfRequiredOutputs(1);
  }
 
}

template< class TGradientImagePixelType, class TTensorPrecision >
typename Image<TGradientImagePixelType,3>::Pointer
DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                    TTensorPrecision >
::GetBaseline()
{
  return static_cast<ScalarImageType* >(this->ProcessObject::GetOutput(1));
}

template< class TGradientImagePixelType, class TTensorPrecision >
void DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                    TTensorPrecision >
::BeforeThreadedGenerateData()
{
  // There need to be at least 6 gradient directions to be able to compute the 
  // tensor basis
  if( m_NumberOfGradientDirections < 7 )
    {
    itkExceptionMacro( << "At least 7 images are required" );
    }
    
  this->ComputeTensorBasis();
}


template< class TGradientImagePixelType, class TTensorPrecision >
void DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                    TTensorPrecision >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int ) 
{
  typename OutputImageType::Pointer outputImage = 
    static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  typename ScalarImageType::Pointer baselineImage = 
    static_cast< ScalarImageType *>(this->ProcessObject::GetOutput(1));
  
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  typedef ImageRegionConstIterator< GradientImagesType > GradientIteratorType;
  typedef typename GradientImagesType::PixelType         GradientVectorType;
  typename GradientImagesType::ConstPointer gradientImagePointer = NULL;
    
  gradientImagePointer =  this->GetInput();
    
  GradientIteratorType git(gradientImagePointer, outputRegionForThread );
  git.GoToBegin();

  typedef ImageRegionIterator< ScalarImageType > ScalarIteratorType;
  ScalarIteratorType bit;
  if(m_EstimateBaseline)
  {
    bit = ScalarIteratorType(baselineImage, outputRegionForThread);
    bit.GoToBegin();
  }

  vnl_vector<TTensorPrecision> B(m_NumberOfGradientDirections);
  vnl_vector<TTensorPrecision> D(7);
    
  while( !git.IsAtEnd() )
  {
    GradientVectorType gv = git.Get();
    for(unsigned int i = 0; i < m_NumberOfGradientDirections; ++i)
      B[i] = gv[i];
    D = EstimateTensor(B);

    TensorPixelType tensor(0.0);
    // First we need to estimate the S_0 then compare it to the threshold
    // D[6] is the estimated S_0
    if( exp(D[6]) >= m_Threshold )
    {                
      // Copy all elements except the estimated S_0 (last element of D)
      std::copy(D.begin(), D.end() - 1, tensor.Begin());

    }

    oit.Set( tensor );
    ++oit; // Output (reconstructed tensor image) iterator
    ++git; // Gradient  image iterator
    if(m_EstimateBaseline)
    {
      bit.Set(static_cast<GradientPixelType>(round(exp(D[6]))));
      ++bit;
    }
    
  }
}


template< class TGradientImagePixelType, class TTensorPrecision >
void DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType, 
                                                    TTensorPrecision >
::ComputeTensorBasis()
{
  if( m_NumberOfGradientDirections < 7 )
    {
    itkExceptionMacro( << "Not enough gradient directions supplied. Need to supply at least 6" );
    }

  m_BMatrix.set_size(m_NumberOfGradientDirections, 7);
  for (unsigned int m = 0; m < m_NumberOfGradientDirections; m++)
    {
    m_BMatrix[m][0] =     -m_BValue * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[0];
    m_BMatrix[m][1] = 2 * -m_BValue * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[1];
    m_BMatrix[m][2] = 2 * -m_BValue * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[2];
    m_BMatrix[m][3] =     -m_BValue * m_GradientDirectionContainer->ElementAt(m)[1] * m_GradientDirectionContainer->ElementAt(m)[1];
    m_BMatrix[m][4] = 2 * -m_BValue * m_GradientDirectionContainer->ElementAt(m)[1] * m_GradientDirectionContainer->ElementAt(m)[2];
    m_BMatrix[m][5] =     -m_BValue * m_GradientDirectionContainer->ElementAt(m)[2] * m_GradientDirectionContainer->ElementAt(m)[2];
    m_BMatrix[m][6] = 1;
    }
 
  m_TensorBasis = vnl_svd<TTensorPrecision>(m_BMatrix).pinverse();
    
}

template< class TGradientImagePixelType, class TTensorPrecision >
void DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType, 
                                                    TTensorPrecision >
::SetGradientImage( GradientDirectionContainerType *gradientDirection, 
                        const GradientImagesType *gradientImage )
{
  this->m_GradientDirectionContainer = gradientDirection;

  m_NumberOfGradientDirections = gradientDirection->Size();

  // ensure that the gradient image we received has as many components as 
  // the number of gradient directions
  if( gradientImage->GetVectorLength() != this->m_NumberOfGradientDirections )
    {
    itkExceptionMacro( << this->m_NumberOfGradientDirections << " gradient directions specified but image has " << gradientImage->GetVectorLength()
      << " components.");
    }
  
  this->ProcessObject::SetNthInput( 0, 
      const_cast< GradientImagesType* >(gradientImage) );
}


template< class TGradientImagePixelType, class TTensorPrecision >
void DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                    TTensorPrecision >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "TensorBasisMatrix: " << m_TensorBasis << std::endl;
  if ( m_GradientDirectionContainer )
    {
    os << indent << "GradientDirectionContainer: "
       << m_GradientDirectionContainer << std::endl;
    }
  else
    {
    os << indent << 
    "GradientDirectionContainer: (Gradient directions not set)" << std::endl;
    }
  os << indent << "NumberOfGradientDirections: " << 
              m_NumberOfGradientDirections << std::endl;
  os << indent << "Threshold for reference B0 image: " << m_Threshold << std::endl;
  os << indent << "BValue: " << m_BValue << std::endl;
}

}

#endif
