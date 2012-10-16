/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNaryVectorAverageImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkNaryVectorAverageImageFilter_h
#define __itkNaryVectorAverageImageFilter_h

#include "itkVectorNaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
  
/** \class NaryVectorAverageImageFilter
 * \brief Implements an operator computing the pixel-wise maximum of several images.
 *
 * This class is parametrized over the types of the input images and the type
 * of the output image.  Numeric conversions (castings) are done by the C++
 * defaults.
 *
 * The pixel type of the output images must have a valid defintion of the
 * operator<. This condition is required because internally this filter will
 * perform an operation similar to:
 *
 *    const OutputPixelType query_value = static_cast<OutputPixelType>(pixel_from_input_n);
 *    if(current_maximum < query_value)
 *      {
 *      current_maximum = query_value; 
 *      }
 * (where current_maximum is also of type OutputPixelType)
 * 
 * for each of the n input images.
 * 
 * For example, this filter could be used directly to find a "maximum projection"
 * of a series of images, often used in preliminary analysis of time-series data.
 *
 * \author Zachary Pincus
 *
 * This filter was contributed by Zachary Pincus from the Department of
 * Biochemistry and Program in Biomedical Informatics at Stanford University
 * School of Medicine
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */

namespace Functor {  
  
template< class TInput, class TOutput >
class VectorAverage
{
public:
  typedef typename NumericTraits< TOutput >::ValueType OutputValueType; 
  // not sure if this typedef really makes things more clear... could just use TOutput?
  
  VectorAverage() {}
  ~VectorAverage() {}
  TOutput operator()( const typename VectorContainer<unsigned int, TInput >::Pointer & B)
    {
      TOutput mean(0.0);

      for(unsigned int i = 0; i < B->size(); ++i)
        {
        mean = mean + B->ElementAt(i);
        }
      mean = mean / B->Size();

      return mean;
    }
  
  bool operator== (const VectorAverage&) const
  {
    return true;
  }
  bool operator!= (const VectorAverage&) const
  {
    return false;
  }
};

}
template <class TInputImage, class TOutputImage>
class ITK_EXPORT NaryVectorAverageImageFilter :
    public
VectorNaryFunctorImageFilter<TInputImage,TOutputImage, 
                             Functor::VectorAverage<  typename TInputImage::PixelType, 
                                                      typename TInputImage::PixelType > > 
{
public:
  /** Standard class typedefs. */
  typedef NaryVectorAverageImageFilter  Self;
  typedef VectorNaryFunctorImageFilter<TInputImage,TOutputImage, 
                                       Functor::VectorAverage< typename TInputImage::PixelType, 
                                                               typename TInputImage::PixelType > >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
protected:
  NaryVectorAverageImageFilter() {}
  virtual ~NaryVectorAverageImageFilter() {}

private:
  NaryVectorAverageImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
