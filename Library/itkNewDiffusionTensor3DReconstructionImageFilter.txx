/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNewDiffusionTensor3DReconstructionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-09-04 20:12:29 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNewDiffusionTensor3DReconstructionImageFilter_txx
#define __itkNewDiffusionTensor3DReconstructionImageFilter_txx

#include "itkNewDiffusionTensor3DReconstructionImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "vnl/vnl_vector.h"

namespace itk {

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
::NewDiffusionTensor3DReconstructionImageFilter()
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 ); 
  m_NumberOfGradientDirections = 0;
  m_Threshold = NumericTraits< ReferencePixelType >::min();
  m_GradientImageTypeEnumeration = Else;
  m_GradientDirectionContainer = NULL;
  m_TensorBasis.set_identity();
  m_BValue = 1.0;
}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
void NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
::BeforeThreadedGenerateData()
{
  // If we have more than 2 inputs, then each input, except the first is a 
  // gradient image. The number of gradient images must match the number of
  // gradient directions.
  const unsigned int numberOfInputs = this->GetNumberOfInputs();

  // There need to be at least 6 gradient directions to be able to compute the 
  // tensor basis
  if( m_NumberOfGradientDirections < 7 )
    {
    itkExceptionMacro( << "At least images are required" );
    }
    
  // If there is only 1 gradient image, it must be an itk::VectorImage. Otherwise 
  // we must have a container of (numberOfInputs-1) itk::Image. Check to make sure
  if ( numberOfInputs == 1
      && m_GradientImageTypeEnumeration != GradientIsInASingleImage )
    {
    std::string gradientImageClassName(
        this->ProcessObject::GetInput(0)->GetNameOfClass());
    if ( strcmp(gradientImageClassName.c_str(),"VectorImage") != 0 )
      {
      itkExceptionMacro( << 
          "There is only one Gradient image. I expect that to be a VectorImage. "
          << "But its of type: " << gradientImageClassName );
      }
    }
    
  this->ComputeTensorBasis();
}


// POTENTIAL WARNING:
//
// Until we fix netlib svd routines, we will need to set the number of thread
// to 1.
template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
void NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int ) 
{
  typename OutputImageType::Pointer outputImage = 
            static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  vnl_vector<double> B(m_NumberOfGradientDirections);
  vnl_vector<double> D(7);
    
  // Two cases here .
  // 1. If the Gradients have been specified in multiple images, we will create
  // 'n' iterators for each of the gradient images and solve the Stejskal-Tanner
  // equations for every pixel. 
  // 2. If the Gradients have been specified in a single multi-component image,
  // one iterator will suffice to do the same.

  if( m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
    std::cerr << "This API no longer functions" << std::endl;
    }
  // The gradients are specified in a single multi-component image
  else if( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
    {
    typedef ImageRegionConstIterator< GradientImagesType > GradientIteratorType;
    typedef typename GradientImagesType::PixelType         GradientVectorType;
    typename GradientImagesType::Pointer gradientImagePointer = NULL;
    
    // Would have liked a dynamic_cast here, but seems SGI doesn't like it
    // The enum will ensure that an inappropriate cast is not done
    gradientImagePointer = static_cast< GradientImagesType * >( 
                              this->ProcessObject::GetInput(0) );
    
    GradientIteratorType git(gradientImagePointer, outputRegionForThread );
    git.GoToBegin();

    while( !git.IsAtEnd() )
      {
      GradientVectorType b = git.Get();
      for(unsigned int i = 0; i < m_NumberOfGradientDirections; ++i)
        {
        if(b[i] == 0)
          B[i] = 0;
        else
          B[i] = log(b[i]);
        }

      TensorPixelType tensor(0.0);

      D = m_TensorBasis * B;

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
      }
    }

}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
void NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
::ComputeTensorBasis()
{
  if( m_NumberOfGradientDirections < 7 )
    {
    itkExceptionMacro( << "Not enough gradient directions supplied. Need to supply at least 6" );
    }

  vnl_matrix<double> bmatrix( m_NumberOfGradientDirections, 7 );
  for (unsigned int m = 0; m < m_NumberOfGradientDirections; m++)
    {
    bmatrix[m][0] =     -m_BValue * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[0];
    bmatrix[m][1] = 2 * -m_BValue * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[1];
    bmatrix[m][2] = 2 * -m_BValue * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[2];
    bmatrix[m][3] =     -m_BValue * m_GradientDirectionContainer->ElementAt(m)[1] * m_GradientDirectionContainer->ElementAt(m)[1];
    bmatrix[m][4] = 2 * -m_BValue * m_GradientDirectionContainer->ElementAt(m)[1] * m_GradientDirectionContainer->ElementAt(m)[2];
    bmatrix[m][5] =     -m_BValue * m_GradientDirectionContainer->ElementAt(m)[2] * m_GradientDirectionContainer->ElementAt(m)[2];
    bmatrix[m][6] = 1;
    }
 
  m_TensorBasis = vnl_svd<double>(bmatrix).pinverse();
    
}

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
void NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
::AddGradientImage( const GradientDirectionType &gradientDirection, 
                    const GradientImageType *gradientImage )
{
  // Make sure crazy users did not call both AddGradientImage and 
  // SetGradientImage
  if( m_GradientImageTypeEnumeration == GradientIsInASingleImage)
    {
    itkExceptionMacro( << "Cannot call both methods:" 
    << "AddGradientImage and SetGradientImage. Please call only one of them.");
    }

  // If the container to hold the gradient directions hasn't been allocated
  // yet, allocate it.
  if( !this->m_GradientDirectionContainer )
    {
    this->m_GradientDirectionContainer = GradientDirectionContainerType::New();
    }
    
  m_GradientDirectionContainer->InsertElement( 
              m_NumberOfGradientDirections, gradientDirection / gradientDirection.two_norm() );
  ++m_NumberOfGradientDirections;
  this->ProcessObject::SetNthInput( m_NumberOfGradientDirections, 
      const_cast< GradientImageType* >(gradientImage) );
  m_GradientImageTypeEnumeration = GradientIsInManyImages;
}

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
void NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
::SetGradientImage( GradientDirectionContainerType *gradientDirection, 
                        const GradientImagesType *gradientImage )
{
  // Make sure crazy users did not call both AddGradientImage and 
  // SetGradientImage
  if( m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
    itkExceptionMacro( << "Cannot call both methods:" 
    << "AddGradientImage and SetGradientImage. Please call only one of them.");
    }

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
  m_GradientImageTypeEnumeration = GradientIsInASingleImage;
}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorPixelType >
void NewDiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorPixelType >
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
  if ( this->m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
    os << indent << "Gradient images haven been supplied " << std::endl;
    }
  else if ( this->m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
    os << indent << "A multicomponent gradient image has been supplied" << std::endl;
    }
}

}

#endif
