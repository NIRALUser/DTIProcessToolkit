/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionWeightedImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-03-26 14:15:33 $
  Version:   $Revision: 1.1 $

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
#include "vnl/algo/vnl_matrix_inverse.h"

namespace itk {

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::DiffusionTensor3DReconstructionWeightedImageFilter()
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 ); 
  m_NumberOfGradientDirections = 0;
  m_NumberOfBaselineImages = 1;
  m_NumberOfIterations = 1;
  m_Threshold = NumericTraits< ReferencePixelType >::min();
  m_GradientImageTypeEnumeration = Else;
  m_GradientDirectionContainer = NULL;
  m_TensorBasis.set_identity();
  m_BValue = 1.0;
}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::BeforeThreadedGenerateData()
{
  // If we have more than 2 inputs, then each input, except the first is a 
  // gradient image. The number of gradient images must match the number of
  // gradient directions.
  const unsigned int numberOfInputs = this->GetNumberOfInputs();

  // There need to be at least 6 gradient directions to be able to compute the 
  // tensor basis
  if( m_NumberOfGradientDirections < 6 )
    {
       itkExceptionMacro( << "At least 6 gradient directions are required" );
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
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int ) 
{

  //PrintSelf(std::cout, indent);

  typename OutputImageType::Pointer outputImage = 
            static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

//  vnl_vector<double> B(m_NumberOfGradientDirections);
//  vnl_vector<double> D(6);
    
  // Two cases here .
  // 1. If the Gradients have been specified in multiple images, we will create
  // 'n' iterators for each of the gradient images and solve the Stejskal-Tanner
  // equations for every pixel. 
  // 2. If the Gradients have been specified in a single multi-component image,
  // one iterator will suffice to do the same.

  
  // The gradients are specified in a single multi-component image
  if( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
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

    // add by Ran    
    typedef ImageRegionConstIterator< TensorImageType > TensorIteratorType;
    //typedef typename TensorImageType::PixelType         TensorVectorType;
    typename TensorImageType::Pointer initialTensorImagePointer = NULL; 
    initialTensorImagePointer = static_cast< TensorImageType * >( 
                              this->ProcessObject::GetInput(1) );
    TensorIteratorType tit(initialTensorImagePointer, outputRegionForThread );
    tit.GoToBegin();


    // Compute the indicies of the baseline images and gradient images
    std::vector<unsigned int> baselineind; // contains the indicies of
                                           // the baseline images
    std::vector<unsigned int> gradientind; // contains the indicies of
                                           // the gradient images

    for(GradientDirectionContainerType::ConstIterator gdcit = this->m_GradientDirectionContainer->Begin();
        gdcit != this->m_GradientDirectionContainer->End(); ++gdcit)
      {
          if(gdcit.Value().one_norm() <= 0.0)
          {
            baselineind.push_back(gdcit.Index());
          }
          else
          {
            gradientind.push_back(gdcit.Index());
          }
      }

    TensorType tensor(0.0);

 //   while( !git.IsAtEnd() )
      while( !git.IsAtEnd() && !tit.IsAtEnd() )
      {
       GradientVectorType b = git.Get();

       typename NumericTraits<ReferencePixelType>::AccumulateType b0 = NumericTraits<ReferencePixelType>::Zero;

       // Average the baseline image pixels
       for(unsigned int i = 0; i < baselineind.size(); ++i)
         {
         b0 += b[baselineind[i]];	//baseline
         }
       b0 /= this->m_NumberOfBaselineImages;
       
       TensorType tensor(tit.Get());
       
       if( (b0 != 0) && (b0 >= m_Threshold) )
         {
         const unsigned int ng = this->m_NumberOfGradientDirections;
           
         vnl_vector<double> s(ng);
         assert(ng == gradientind.size());
         for(unsigned int i = 0; i < ng; ++i)
           {
           if(b[gradientind[i]])
             s[i] = log(b[gradientind[i]]) - log(b0);
           else
             s[i] = - log(b0);
           }
         
//         std::cout << "Orig tensor: " << tensor << std::endl;
         for(unsigned int iter = 0; iter < m_NumberOfIterations; ++iter)
           {
           vnl_vector<double> phi(ng), tmp(ng), beta(6);
           std::copy(tensor.Begin(),tensor.End(),beta.begin());
           tmp = -m_BValue * m_BMatrix * beta;
           
           for(unsigned int i = 0; i < ng; ++i)
             {
             phi[i] = exp(tmp[i]);
             phi[i] *= phi[i];
             }
         
           vnl_diag_matrix<double> W2(phi);
           vnl_vector<double> wb = 
             vnl_svd<double>(-m_BValue*m_BMatrix.transpose() * W2 * -m_BValue * m_BMatrix).solve(-m_BValue * m_BMatrix.transpose()*W2*s);
           
           std::copy(wb.begin(),wb.end(),tensor.Begin());

//           std::cout << "Iterations [" << iter << "]: " << tensor <<std::endl;
           }

         }
      
       oit.Set( tensor );
       // oit.Set( tensor / m_BValue );
       ++oit; // Output (reconstructed tensor image) iterator
       ++git; // Gradient  image iterator
       ++tit; // initial tensor, add by Ran
      }
    }

}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::ComputeTensorBasis()
{
  if( m_NumberOfGradientDirections < 6 )
    {
    itkExceptionMacro( << "Not enough gradient directions supplied. Need to supply at least 6" );
    }

  // This is only important if we are using a vector image.  For
  // images added one at a time, this is not needed but doesn't hurt.
  std::vector<unsigned int> gradientind;
  for(GradientDirectionContainerType::ConstIterator gdcit = this->m_GradientDirectionContainer->Begin();
      gdcit != this->m_GradientDirectionContainer->End(); ++gdcit)
    {
    if(gdcit.Value().one_norm() > 0.0)
      {
        gradientind.push_back(gdcit.Index());
      }
    }

  m_BMatrix.set_size( m_NumberOfGradientDirections, 6 );
  for (unsigned int m = 0; m < m_NumberOfGradientDirections; m++)
    {
    m_BMatrix[m][0] =     m_GradientDirectionContainer->ElementAt(gradientind[m])[0] * m_GradientDirectionContainer->ElementAt(gradientind[m])[0];
    m_BMatrix[m][1] = 2 * m_GradientDirectionContainer->ElementAt(gradientind[m])[0] * m_GradientDirectionContainer->ElementAt(gradientind[m])[1];
    m_BMatrix[m][2] = 2 * m_GradientDirectionContainer->ElementAt(gradientind[m])[0] * m_GradientDirectionContainer->ElementAt(gradientind[m])[2];
    m_BMatrix[m][3] =     m_GradientDirectionContainer->ElementAt(gradientind[m])[1] * m_GradientDirectionContainer->ElementAt(gradientind[m])[1];
    m_BMatrix[m][4] = 2 * m_GradientDirectionContainer->ElementAt(gradientind[m])[1] * m_GradientDirectionContainer->ElementAt(gradientind[m])[2];
    m_BMatrix[m][5] =     m_GradientDirectionContainer->ElementAt(gradientind[m])[2] * m_GradientDirectionContainer->ElementAt(gradientind[m])[2];
    }

  m_TensorBasis = vnl_svd<double>(m_BMatrix).pinverse();
}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::SetGradientImage( GradientDirectionContainerType *gradientDirection, 
                         GradientImagesType *gradientImage )
{
  // Make sure crazy users did not call both AddGradientImage and 
  // SetGradientImage
  if( m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
    itkExceptionMacro( << "Cannot call both methods:" 
    << "AddGradientImage and SetGradientImage. Please call only one of them.");
    }

  this->m_GradientDirectionContainer = gradientDirection;

  unsigned int numImages = gradientDirection->Size();
  this->m_NumberOfBaselineImages = 0;
  for(GradientDirectionContainerType::Iterator it = this->m_GradientDirectionContainer->Begin();
      it != this->m_GradientDirectionContainer->End(); it++)
    {
    if(it.Value().one_norm() <= 0.0)
      {
        this->m_NumberOfBaselineImages++;
      }
    else // Normalize non-zero gradient directions
      {
        it.Value() = it.Value() / it.Value().two_norm();
      }
    }
      
  this->m_NumberOfGradientDirections = numImages - this->m_NumberOfBaselineImages;

  // ensure that the gradient image we received has as many components as 
  // the number of gradient directions
  if( gradientImage->GetVectorLength() != this->m_NumberOfBaselineImages + this->m_NumberOfGradientDirections )
    {
    itkExceptionMacro( << this->m_NumberOfGradientDirections << " gradients + " << this->m_NumberOfBaselineImages
                       << "baselines = " << this->m_NumberOfGradientDirections + this->m_NumberOfBaselineImages
                       << " directions specified but image has " << gradientImage->GetVectorLength()
      << " components.");
    }
  
  this->ProcessObject::SetNthInput( 0, 
      const_cast< GradientImagesType* >(gradientImage) );
  m_GradientImageTypeEnumeration = GradientIsInASingleImage;
}

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::SetInitialTensor( TensorImageType *tensorImage )
{
  this->ProcessObject::SetNthInput( 1, 
      const_cast< TensorImageType* >(tensorImage) );
}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionWeightedImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "TensorBasisMatrix: " << m_TensorBasis << std::endl;
  os << indent << "Coeffs: " << m_BMatrix << std::endl;
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
  os << indent << "NumberOfBaselineImages: " << 
              m_NumberOfBaselineImages << std::endl;
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
