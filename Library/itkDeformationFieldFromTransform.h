/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldFromTransform.h,v $
  Language:  C++
  Date:      $Date: 2007-09-04 20:12:29 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDeformationFieldFromTransform_h
#define __itkDeformationFieldFromTransform_h

#include "itkImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTransform.h"

namespace itk
{

/** \class DeformationFieldFromTransform
 * \brief Computes a deformation field from two sets of landmarks.
 *
 * DeformationFieldFromTransform produces a deformation field from two set of input
 * landmarks.  One set of landmarks are associated to the input space while the
 * second set of landmarks is associated with the output space.
 *
 * A KernelBase spline is used to interpolate the deformations and produce
 * deformation values for all the nodes of the image grid that will be produced
 * as output.
 * 
 * The number of landmarks in the KernelBased spline will have a dramatic
 * effect on both the precision of output deformation field and the
 * computational time required for the filter to complete the estimation. 
 *
 * 
 * This source object expects the image to be of pixel type Vector.
 *
 * \ingroup ImageSource
 */
template <class TOutputImage, class TPrecision>
class ITK_EXPORT DeformationFieldFromTransform:
    public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DeformationFieldFromTransform         Self;
  typedef ImageSource<TOutputImage>      Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;
  
  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(DeformationFieldFromTransform, ImageSource);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Image size typedef. */
  typedef typename OutputImageType::SizeType      OutputSizeType;

  /** Image index typedef. */
  typedef typename OutputImageType::IndexType     OutputIndexType;

  /** Image pixel value typedef. */
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename OutputPixelType::ValueType     OutputPixelComponentType;

  /** Image spacing typedef */
  typedef typename TOutputImage::SpacingType      SpacingType;
  typedef typename TOutputImage::PointType        OriginPointType;

  /** Precision used for transform */
  typedef TPrecision PrecisionType;
  
  /** Transform typedef */
  typedef Transform<PrecisionType,
                    TOutputImage::ImageDimension,
                    TOutputImage::ImageDimension> TransformType;
  
  /** Set the transform. */
  itkSetObjectMacro( Transform, TransformType );

  /** Get the transform. */
  itkGetObjectMacro( Transform, TransformType );
  
  /** Set the size of the output image. */
  itkSetMacro( OutputRegion, OutputImageRegionType );

  /** Get the size of the output image. */
  itkGetConstReferenceMacro( OutputRegion, OutputImageRegionType );
     
  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void SetOutputSpacing( const double* values);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, OriginPointType);
  virtual void SetOutputOrigin( const double* values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro( OutputOrigin, OriginPointType );

  /** DeformationFieldFromTransform produces an image. As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation();


  /** Method Compute the Modified Time based on changed to the components. */
  unsigned long GetMTime( void ) const;

protected:
  DeformationFieldFromTransform();
  ~DeformationFieldFromTransform() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** 
   * GenerateData() computes the internal KernelBase spline and resamples
   * the deformation field.
   */
  void GenerateData();

private:
  DeformationFieldFromTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  OutputImageRegionType         m_OutputRegion;      // Region of the output image
  SpacingType                   m_OutputSpacing;     // output image spacing
  OriginPointType               m_OutputOrigin;      // output image origin
  
  typename TransformType::Pointer  m_Transform;         // transform object
};

  
} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformationFieldFromTransform.txx"
#endif
  
#endif
