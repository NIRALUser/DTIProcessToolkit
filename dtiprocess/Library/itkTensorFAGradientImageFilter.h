/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorFAGradientImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorFAGradientImageFilter_h
#define __itkTensorFAGradientImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkDiffusionTensor3D.h>

namespace itk
{


/** \class TensorFAGradientImageFilter
 * \brief Computes the 6-element vector field that are the unique
 * elements of the matrix-logarithm of the tensor field.
 *
 * ExpEuclideanImageFilter applies pixel-wise the invokation for
 * computing the matrix logarithm of every pixel. 
 * 
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template<typename T>
class ITK_EXPORT TensorFAGradientImageFilter :
    public
ImageToImageFilter<Image<DiffusionTensor3D<T>, 3>,
                   Image<CovariantVector<T,3>, 3> >
{
public:
//  typedef Image<Vector<typename TInputImage::PixelType::RealValueType,6>,3 > TOutputImage;
  typedef Image<DiffusionTensor3D<T>, 3> InputImageType;
  typedef Image<CovariantVector<T,3>, 3> OutputImageType;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef T RealType;

  /** Standard class typedefs. */
  typedef TensorFAGradientImageFilter  Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  itkSetMacro(Sigma, double);

  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  
  virtual void GenerateData();


protected:
  TensorFAGradientImageFilter() {};
  virtual ~TensorFAGradientImageFilter() {};

private:
  TensorFAGradientImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double m_Sigma;

};


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorFAGradientImageFilter.txx"
#endif

}


#endif
