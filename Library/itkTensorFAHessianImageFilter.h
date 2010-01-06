/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorFAHessianImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorFAHessianImageFilter_h
#define __itkTensorFAHessianImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkDiffusionTensor3D.h>

#include "tensorderivs.h"

namespace itk
{


/** \class TensorFAHessianImageFilter
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
class ITK_EXPORT TensorFAHessianImageFilter :
    public
ImageToImageFilter<Image<DiffusionTensor3D<T>, 3>,
                   Image<SymmetricSecondRankTensor<3,T>, 3> >
{
public:
//  typedef Image<Vector<typename TInputImage::PixelType::RealValueType,6>,3 > TOutputImage;
  typedef Image<DiffusionTensor3D<T>, 3> InputImageType;
  typedef Image<SymmetricSecondRankTensor<3,T>, 3> InputImageType;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename OutputPixelType::EigenVectorsMatrixType EigenVectorType;
  typedef typename OutputPixelType::EigenValuesArrayType EigenValueType;
  typedef T RealType;


  /** Standard class typedefs. */
  typedef TensorFAHessianImageFilter  Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  
  virtual void GenerateData();


protected:
  TensorFAHessianImageFilter() {};
  virtual ~TensorFAHessianImageFilter() {};

private:
  TensorFAHessianImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double m_Scale;

};


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorFAHessianImageFilter.txx"
#endif

}


#endif
