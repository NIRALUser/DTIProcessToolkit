/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLinearInterpolateTensorImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLinearInterpolateTensorImageFunction_h
#define _itkLinearInterpolateTensorImageFunction_h

#include "itkInterpolateImageFunction.h"

namespace itk
{

/** \class LinearInterpolateTensorImageFunction
 * \brief Linearly interpolate an image at specified positions.
 *
 * LinearInterpolateTensorImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type 
 * (e.g. float or double).
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for images with scalar pixel
 * types. For vector images use VectorLinearInterpolateTensorImageFunction.
 *
 * \sa VectorLinearInterpolateTensorImageFunction
 *
 * \ingroup ImageFunctions ImageInterpolators 
 */
template <class TInputImage, class TCoordRep = float>
class ITK_EXPORT LinearInterpolateTensorImageFunction : 
  public InterpolateTensorImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef LinearInterpolateTensorImageFunction Self;
  typedef InterpolateTensorImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(LinearInterpolateTensorImageFunction, InterpolateTensorImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** PixelType typedef supprt. */
  typedef typename InputImageType::PixelType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const;

protected:
  LinearInterpolateTensorImageFunction();
  ~LinearInterpolateTensorImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  LinearInterpolateTensorImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  

};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLinearInterpolateTensorImageFunction.txx"
#endif

#endif
