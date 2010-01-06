/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorDotProductImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorDotProductImageFilter_h
#define __itkVectorDotProductImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"


namespace itk
{
  
/** \class VectorDotProductImageFilter
 * \brief Implements an operator for pixel-wise masking of the input 
 * image with the negative of a mask.
 *
 * This class is parametrized over the types of the  
 * input image type, the mask image type and the type of the output image. 
 * Numeric conversions (castings) are done by the C++ defaults.
 *
 * The pixel type of the input 2 image must have a valid defintion of the
 * operator != with zero . This condition is required because internally this
 * filter will perform the operation
 *
 *        if pixel_from_mask_image != 0 
 *             pixel_output_image = 0
 *        else
 *             pixel_output_image = pixel_input_image
 *
 * The pixel from the input 1 is cast to the pixel type of the output image.
 *
 * Note that the input and the mask images must be of the same size.
 *
 * \warning Any pixel value other than 0 will not be masked out. 
 *
 * \sa DotProductImageFilter
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  
  
template< class TInput1, class TInput2, class TOutput >
class VectorDotProductInput
{
public:
  typedef typename NumericTraits< TInput1 >::AccumulateType AccumulatorType;

  VectorDotProductInput() {};
  ~VectorDotProductInput() {};
  bool operator!=( const VectorDotProductInput & ) const
  {
    return false;
  }
  bool operator==( const VectorDotProductInput & other ) const
  {
    return !(*this != other);
  }
  TOutput operator()( const TInput1 & A, const TInput2 & B)
  {
    if (B == NumericTraits< TInput2 >::Zero ) 
      {
      TInput1 R = A;
      R.Fill(0);
      return R;
      }
    else
      {
      return static_cast<TOutput>( A );
      }
  }
}; 

}
template <class TInput1Image, class TInput2Image, class TOutputImage>
class ITK_EXPORT VectorDotProductImageFilter :
    public
BinaryFunctorImageFilter<TInput1Image,TInput2Image,TOutputImage, 
                         Functor::VectorDotProductInput< 
  typename TInput1Image::PixelType, 
  typename TInput2Image::PixelType,
  typename TOutputImage::PixelType>   >


{
public:
  /** Standard class typedefs. */
  typedef VectorDotProductImageFilter  Self;
  typedef BinaryFunctorImageFilter<TInput1Image,TInput2Image,TOutputImage, 
                                   Functor::VectorDotProductInput< 
    typename TInput1Image::PixelType, 
    typename TInput2Image::PixelType,
    typename TOutputImage::PixelType>   
  >  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
protected:
  VectorDotProductImageFilter() {}
  virtual ~VectorDotProductImageFilter() {}

private:
  VectorDotProductImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
