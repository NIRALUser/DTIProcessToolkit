/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionRicianImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionRicianImageFilter_txx
#define __itkDiffusionTensor3DReconstructionRicianImageFilter_txx

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "vnl/vnl_vector.h"
#include "itkDiffusionTensor3DReconstructionRicianImageFilter.h"
#include "vnl/algo/vnl_lbfgs.h"
#include "vnl/algo/vnl_powell.h"
#include "vnl/algo/vnl_conjugate_gradient.h"
#include "itkLBFGSBOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
//#include "itkRegularStepGradientDescentOptimizer.h"

namespace itk {

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::DiffusionTensor3DReconstructionRicianImageFilter()
{
  // At least 1 inputs is necessary for a vector image.
  // For images added one at a time we need at least six
  this->SetNumberOfRequiredInputs( 1 ); 
  m_NumberOfGradientDirections = 0;
  m_NumberOfBaselineImages = 1;
  m_Threshold = NumericTraits< ReferencePixelType >::min();
  m_GradientImageTypeEnumeration = Else;
  m_GradientDirectionContainer = NULL;
  m_TensorBasis.set_identity();
  m_BValue = 1.0;
}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
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
void DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int ) 
{
  typename OutputImageType::Pointer outputImage = 
            static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  vnl_vector<double> B(m_NumberOfGradientDirections);
  vnl_vector<double> D(6);
    
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

//    while( !git.IsAtEnd() )
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

     // TensorType tensor(0.0);
      tensor = tit.Get();

      if( (b0 != 0) && (b0 >= m_Threshold) )
        {

         const unsigned int ng = this->m_NumberOfGradientDirections;
         vnl_vector<ReferencePixelType> s(ng);

//         std::cout << "b: " << b << std::endl;
         std::cout << "ind: " << tit.GetIndex() << std::endl;
        //B is the value corresponding to the gradient direction
         for( unsigned int i = 0; i< m_NumberOfGradientDirections; i++ )
           s[i] = b[gradientind[i]];
         
          
         typedef RicianLikelihood<ReferencePixelType> FittingFunctionType;

         typename FittingFunctionType::Pointer loglikelihood = FittingFunctionType::New();
         
         loglikelihood->SetS0(b0);
         loglikelihood->SetDesign_Matrix(-m_BValue * m_BMatrix);
         loglikelihood->SetSigma(m_Sigma);
         loglikelihood->SetSignal(s);

         Array<long int> bounds(6);
         bounds.Fill(0);

         Array<double> lbounds(6);
         lbounds.Fill(-1);
//         lbounds[0] = lbounds[3] = lbounds[5] = 0.0;
         
         Array<double> hbounds(6);
         hbounds.Fill(1);

         Array<double> vnlt(tensor.GetDataPointer(),6);

//         std::cout << "init: " << vnlt << std::endl;
//         typedef LBFGSBOptimizer OptimizerType;
         typedef RegularStepGradientDescentOptimizer OptimizerType;
         OptimizerType::Pointer optimizer = OptimizerType::New();
         optimizer->SetCostFunction(loglikelihood);
         optimizer->SetInitialPosition(vnlt);
//          optimizer->SetBoundSelection(bounds);
//          optimizer->SetLowerBound(lbounds);
//          optimizer->SetUpperBound(hbounds);
         optimizer->MinimizeOn();
//         optimizer->SetCostFunctionConvergenceFactor(1.0e7);
         optimizer->SetMaximumStepLength(1.0e-4);
         optimizer->SetMinimumStepLength(1.0e-15);
         optimizer->SetRelaxationFactor(.2);

//         optimizer->SetLearningRate(1.0e-8);
//         optimizer->DebugOn();
//         std::cout << "Initial Cost: "<< loglikelihood->GetValue(vnlt) << std::endl;

         try
           {
           optimizer->StartOptimization();        
           vnlt = optimizer->GetCurrentPosition();
           } 
         catch(itk::ExceptionObject &e)
           {
           std::cerr << e << std::endl;
           }

//         std::cout << "min: " << vnlt << std::endl;
//         std::cout << "Final Cost: "<< optimizer->GetValue() << std::endl;
           
         std::copy(vnlt.begin(),vnlt.end(),tensor.Begin());

        }

      oit.Set( tensor );
  //  oit.Set( tensor / m_BValue);
      ++oit; // Output (reconstructed tensor image) iterator
      ++git; // Gradient  image iterator
      ++tit; // initial tensor, add by Ran
      }
    }

}


template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
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
void DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
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
void DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
  TGradientImagePixelType, TTensorType >
::SetInitialTensor( TensorImageType *tensorImage )
{
  this->ProcessObject::SetNthInput( 1, 
      const_cast< TensorImageType* >(tensorImage) );
}

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType, class TTensorType >
void DiffusionTensor3DReconstructionRicianImageFilter< TReferenceImagePixelType,
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

} // end namespace itk
#endif
