/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorNaryFunctorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorNaryFunctorImageFilter_h
#define __itkVectorNaryFunctorImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkImageIterator.h"
#include "itkVectorContainer.h"

namespace itk
{
  
/** \class VectorNaryFunctorImageFilter
 * \brief Implements pixel-wise generic operation of Nth similar images.
 *
 * This class is parameterized over the types of the input images
 * and the type of the output image.  It is also parameterized by the
 * operation to be applied.  A Functor style is used to represent the
 * function.
 *
 * All the input images are of the same type.
 * 
 * \ingroup IntensityImageFilters   Multithreaded
 */

template <class TInputImage, class TOutputImage, class TFunction >
class ITK_EXPORT VectorNaryFunctorImageFilter :
    public InPlaceImageFilter<TInputImage,TOutputImage> 

{
public:
  /** Standard class typedefs. */
  typedef VectorNaryFunctorImageFilter  Self;
  typedef InPlaceImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorNaryFunctorImageFilter, InPlaceImageFilter);

  /** Some typedefs. */
  typedef TFunction   FunctorType;
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType; 
  typedef typename InputImageType::PixelType    InputImagePixelType; 
  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;
  typedef typename OutputImageType::PixelType   OutputImagePixelType;
  typedef VectorContainer<unsigned int, InputImagePixelType > VectorNaryArrayType; 

  /** Get the functor object.  The functor is returned by reference.
   * (Functors do not have to derive from itk::LightObject, so they do
   * not necessarily have a reference count. So we cannot return a
   * SmartPointer). */
  FunctorType& GetFunctor() { return m_Functor; };

  /** Set the functor object.  This replaces the current Functor with a
   * copy of the specified Functor. This allows the user to specify a
   * functor that has ivars set differently than the default functor.
   * This method requires an operator!=() be defined on the functor
   * (or the compiler's default implementation of operator!=() being
   * appropriate). */
  void SetFunctor(FunctorType& functor)
  {
    if ( m_Functor != functor )
      {
      m_Functor = functor;
      this->Modified();
      }
  }
  
protected:
  VectorNaryFunctorImageFilter();
  virtual ~VectorNaryFunctorImageFilter() {};

  /** VectorNaryFunctorImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData() routine
   * which is called for each processing thread. The output image data is
   * allocated automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to the
   * portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            int threadId );

private:
  VectorNaryFunctorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  FunctorType m_Functor;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorNaryFunctorImageFilter.txx"
#endif

#endif
