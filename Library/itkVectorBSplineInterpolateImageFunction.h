/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorBSplineInterpolateImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorBSplineInterpolateImageFunction_h
#define __itkVectorBSplineInterpolateImageFunction_h

#include "itkVectorInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace itk
{

/** 
 * \class VectorBSplineInterpolateImageFunction
 * \brief BSplinely interpolate a vector image at specified positions.
 *
 * VectorBSplineInterpolateImageFunction linearly interpolates a vector
 * image intensity non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type.
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for Vector images. For
 * scalar images use BSplineInterpolateImageFunction.
 *
 * \ingroup ImageFunctions ImageInterpolators
 * 
 */
template <class TInputImage, class TCoordRep = float, class TCoefficientType = double>
class ITK_EXPORT VectorBSplineInterpolateImageFunction : 
  public VectorInterpolateImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef VectorBSplineInterpolateImageFunction Self;
  typedef VectorInterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorBSplineInterpolateImageFunction, 
    VectorInterpolateImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::ValueType      ValueType;
  typedef typename Superclass::RealType       RealType;
    
  /** Grab the vector dimension from the superclass. */
  itkStaticConstMacro(Dimension, unsigned int,
                       Superclass::Dimension);

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is Vector<double,Dimension> */
  typedef typename Superclass::OutputType OutputType;

  /** component extractor type */
  typedef typename PixelType::ComponentType  ComponentType;
  typedef Image<ComponentType, ImageDimension> ComponentImageType;
  typedef VectorIndexSelectionCastImageFilter<InputImageType,ComponentImageType> ComponentAdaptorType;
  typedef typename ComponentAdaptorType::Pointer ComponentAdaptorPointer;

  /** component interpolator type */
  typedef BSplineInterpolateImageFunction<ComponentImageType, TCoordRep, TCoefficientType> ComponentInterpolateFunctionType;
  typedef typename ComponentInterpolateFunctionType::Pointer ComponentInterpolateFunctionPointer;

  /** Set the input image.  This must be set by the user. */
  virtual void SetInputImage(const TInputImage * inputData);

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the spline interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const;

protected:
  VectorBSplineInterpolateImageFunction();
  ~VectorBSplineInterpolateImageFunction(){};
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:
  VectorBSplineInterpolateImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  
 
  std::vector<ComponentAdaptorPointer> m_ComponentAdaptors;
  std::vector<ComponentInterpolateFunctionPointer> m_ComponentInterpolators;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorBSplineInterpolateImageFunction.txx"
#endif

#endif
