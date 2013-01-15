/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldFromTransform.txx,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#ifndef _itkDeformationFieldFromTransform_txx
#define _itkDeformationFieldFromTransform_txx

#include "itkDeformationFieldFromTransform.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{

/**
 * Initialize new instance
 */
template <class TOutputImage, class TPrecision>
DeformationFieldFromTransform<TOutputImage, TPrecision>
::DeformationFieldFromTransform()
{
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);
  
}


/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template <class TOutputImage, class TPrecision>
void 
DeformationFieldFromTransform<TOutputImage, TPrecision>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "OutputRegion:    " << m_OutputRegion << std::endl;
  os << indent << "OutputSpacing:   " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin:    " << m_OutputOrigin << std::endl;
}



/**
 * Set the output image spacing.
 */
template <class TOutputImage, class TPrecision>
void 
DeformationFieldFromTransform<TOutputImage, TPrecision>
::SetOutputSpacing(
  const double* spacing)
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}


/**
 * Set the output image origin.
 */
template <class TOutputImage, class TPrecision>
void 
DeformationFieldFromTransform<TOutputImage, TPrecision>
::SetOutputOrigin(
  const double* origin)
{
  OriginPointType p(origin);
  this->SetOutputOrigin( p );
}


/**
 * ThreadedGenerateData
 */
template <class TOutputImage, class TPrecision>
void 
DeformationFieldFromTransform<TOutputImage, TPrecision>
::ThreadedGenerateData(const OutputImageRegionType& outputRegion,
                       int threadId)
{

  itkDebugMacro(<<"Actually executing");

  // Get the output pointers
  OutputImageType *  outputPtr = this->GetOutput();

  // Create an iterator that will walk the output region for this thread.
  typedef ImageRegionIteratorWithIndex< 
                                  TOutputImage> OutputIterator;

  OutputIterator outIt( outputPtr, outputRegion );

  // Define a few indices that will be used to translate from an input pixel
  // to an output pixel
  OutputIndexType outputIndex;         // Index to current output pixel

  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;

  InputPointType outputPoint;    // Coordinates of current output pixel

  // Support for progress methods/callbacks
  ProgressReporter progress(this, 0, outputRegion.GetNumberOfPixels(), 10);
        
  outIt.GoToBegin();

  // Walk the output region
  while ( !outIt.IsAtEnd() )
    {
    // Determine the index of the current output pixel
    outputIndex = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint( outputIndex, outputPoint );

    // Compute corresponding inverse displacement vector
    OutputPointType interpolatedDeformation = 
                        m_Transform->TransformPoint( outputPoint );

    OutputPixelType displacement;
    for( unsigned int i=0; i < ImageDimension; i++)
      {
      displacement[i] = interpolatedDeformation[i] - outputPoint[i];

      }

    outIt.Set( displacement );
    ++outIt;
    progress.CompletedPixel();
    }

  return;
}


/** 
 * Inform pipeline of required output region
 */
template <class TOutputImage, class TPrecision>
void 
DeformationFieldFromTransform<TOutputImage, TPrecision>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  OutputImagePointer outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }

  // Set the size of the output region
  outputPtr->SetLargestPossibleRegion( m_OutputRegion );

  // Set spacing and origin
  outputPtr->SetSpacing( m_OutputSpacing );
  outputPtr->SetOrigin( m_OutputOrigin );

  return;
}



/** 
 * Verify if any of the components has been modified.
 */
template <class TOutputImage, class TPrecision>
unsigned long 
DeformationFieldFromTransform<TOutputImage, TPrecision>
::GetMTime( void ) const
{
  unsigned long latestTime = Object::GetMTime(); 

  return latestTime;
}



} // end namespace itk

#endif
