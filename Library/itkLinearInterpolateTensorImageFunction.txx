/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLinearInterpolateTensorImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLinearInterpolateTensorImageFunction_txx
#define _itkLinearInterpolateTensorImageFunction_txx
#include "itkLinearInterpolateTensorImageFunction.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep>
const unsigned long
LinearInterpolateTensorImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep>
LinearInterpolateTensorImageFunction< TInputImage, TCoordRep >
::LinearInterpolateTensorImageFunction()
{

}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void
LinearInterpolateTensorImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename LinearInterpolateTensorImageFunction< TInputImage, TCoordRep >
::OutputType
LinearInterpolateTensorImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType& index) const
{
  unsigned int dim;  // index over dimension

  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  signed long baseIndex[ImageDimension];
  double distance[ImageDimension];
  long tIndex;

  for( dim = 0; dim < ImageDimension; dim++ )
    {
    // The following "if" block is equivalent to the following line without
    // having to call floor.
    //    baseIndex[dim] = (long) floor( index[dim] );
    if (index[dim] >= 0.0)
      {
      baseIndex[dim] = (long) index[dim];
      }
    else
      {
      tIndex = (long) index[dim];
      if (double(tIndex) != index[dim])
        {
        tIndex--;
        }
      baseIndex[dim] = tIndex;
      }
    distance[dim] = index[dim] - double( baseIndex[dim] );
    }
  
  /**
   * Interpolated value is the weighted sum of each of the surrounding
   * neighbors. The weight for each neighbor is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  itk::VectorContainer<unsigned int, double> weights;
  itk::VectorContainer<unsigned int, PixelType> pixels;

  for( unsigned int counter = 0; counter < m_Neighbors; counter++ )
    {

    double overlap = 1.0;          // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour
    IndexType neighIndex;

    // get neighbor index and overlap fraction
    for( dim = 0; dim < ImageDimension; dim++ )
      {

      if ( upper & 1 )
        {
        neighIndex[dim] = baseIndex[dim] + 1;
        overlap *= distance[dim];
        }
      else
        {
        neighIndex[dim] = baseIndex[dim];
        overlap *= 1.0 - distance[dim];
        }

      upper >>= 1;

      }
    
    // get neighbor value only if overlap is not zero
    if( overlap )
      {
      weights.push_back(overlap);
      pixels.push_back(this->GetInputImage()->GetPixel(neighIndex));
      totalOverlap += overlap;
      }

    if( totalOverlap == 1.0 )
      {
      // finished
      break;
      }

    }

  PixelType mean;
  SymmetricSpaceGeometry<double> ssg;
  TensorStatistics<double> ts(ssg);
  
  ts.ComputeWeightedAve(weights,pixels,mean);

  return mean;
}

} // namespace itk

#endif
