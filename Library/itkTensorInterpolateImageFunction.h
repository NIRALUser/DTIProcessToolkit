/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorInterpolateImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorInterpolateImageFunction_h
#define __itkTensorInterpolateImageFunction_h

#include "itkImageFunction.h"
#include "itkFixedArray.h"

namespace itk
{

/**
 * Due to a bug in MSVC, an enum value cannot be accessed out of a template
 * parameter until the template class opens.  In order for templated classes
 * to access the dimension of a template parameter in defining their
 * own dimension, this class is needed as a work-around.
 */
// template <typename T>
// struct GetDimension
// {
//   itkStaticConstMacro(Dimension, int, T::Dimension);
// }; 

  
/** \class TensorInterpolateImageFunction
 * \brief Base class for all vector image interpolaters.
 *
 * TensorInterpolateImageFunction is the base for all ImageFunctions that
 * interpolates image with vector pixel types. This function outputs
 * a return value of type Tensor<double,Dimension>.
 *
 * This class is templated input image type and the coordinate
 * representation type.
 *
 * \warning This hierarchy of functions work only for images 
 * with Tensor-based pixel types. For scalar images use 
 * InterpolateImageFunction.
 * 
 * \sa InterpolateImageFunction
 * \ingroup ImageFunctions ImageInterpolators
 */
template <
class TInputImage,
class TCoordRep = float,
class TPixelType = ITK_TYPENAME TInputImage::PixelType
>
class ITK_EXPORT TensorInterpolateImageFunction : 
  public ImageFunction<
    TInputImage, 
    TPixelType,
    TCoordRep > 
{
public:
  /** Extract the vector dimension from the pixel template parameter. */
  itkStaticConstMacro(Dimension, unsigned int,
                      TPixelType::InternalDimension);
  
  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef TensorInterpolateImageFunction Self;
  typedef ImageFunction<TInputImage,
    TPixelType, TCoordRep > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorInterpolateImageFunction, ImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename InputImageType::PixelType  PixelType;
  typedef typename PixelType::ValueType       ValueType;
  typedef typename NumericTraits<ValueType>::RealType  RealType;
    

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is FixedArray<RealType,Dimension>. */
  typedef typename Superclass::OutputType OutputType;

  /** CoordRep typedef support. */
  typedef TCoordRep CoordRepType;

  /** Returns the interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType Evaluate( const PointType& point ) const
    {
    ContinuousIndexType index;
    this->GetInputImage()->TransformPhysicalPointToContinuousIndex( point, index );
    return ( this->EvaluateAtContinuousIndex( index ) );
    }

  /** Interpolate the image at a continuous index position
   *
   * Returns the interpolated image intensity at a 
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * Subclasses must override this method.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const = 0;

  /** Interpolate the image at an index position.
   * Simply returns the image value at the
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const
    {
    OutputType output;
    PixelType input = this->GetInputImage()->GetPixel( index );
    for( unsigned int k = 0; k < Dimension; k++ )
      {
      output[k] = static_cast<double>( input[k] );
      }
    return ( output );
    }

protected:
  TensorInterpolateImageFunction() {}
  ~TensorInterpolateImageFunction() {}
  void PrintSelf(std::ostream& os, Indent indent) const
    { Superclass::PrintSelf( os, indent ); }

private:
  TensorInterpolateImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif


