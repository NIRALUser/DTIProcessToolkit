/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkExpEuclideanTensorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008-04-11 16:31:05 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkExpEuclideanTensorImageFilter_h
#define __itkExpEuclideanTensorImageFilter_h

#include <itkImageToImageFilter.h>
#include <vnl/vnl_double_3x3.h>

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
// namespace Functor {  
 
// template< typename TInput >
// class ExpEuclideanTensorFunction
// {
// public:
//   typedef typename TInput::RealValueType  RealValueType;
//   typedef typename itk::Vector<RealValueType, 6> OutputType;
//   typedef TInput TensorType;
//   typedef typename TInput::EigenVectorsMatrixType EigenVectorType;
//   typedef typename TInput::EigenValuesArrayType EigenValueType;

//   ExpEuclideanTensorFunction() {}
//   ~ExpEuclideanTensorFunction() {}
//   inline OutputType operator()( const TInput & x )
//   {
//     EigenValueType D;
//     EigenVectorType U;

//     x.ComputeEigenAnalysis(D,U);
    
//     itk::Matrix<RealValueType,3,3> diag;
    
//     diag(0,0) = log(D[0]);
//     diag(1,1) = log(D[1]);
//     diag(2,2) = log(D[2]);

//     itk::Matrix<RealValueType,3,3> matlog;
//     matlog = U * diag * U.GetTranspose();

//     OutputType result;
//     result[0] = matlog(0,0);
//     result[1] = matlog(0,1);
//     result[2] = matlog(0,2);
//     result[3] = matlog(1,1);
//     result[4] = matlog(1,2);
//     result[5] = matlog(2,2);

//     return result;
//   }
// }; 


// }  // end namespace functor


/** \class ExpEuclideanTensorImageFilter
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
class ITK_EXPORT ExpEuclideanTensorImageFilter :
    public
ImageToImageFilter<Image<Vector<T,6>, 3>,
                   Image<DiffusionTensor3D<T>, 3> >
{
public:
//  typedef Image<Vector<typename TInputImage::PixelType::RealValueType,6>,3 > TOutputImage;
  typedef Image<Vector<T,6>, 3> InputImageType;
  typedef Image<DiffusionTensor3D<T>, 3> OutputImageType;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename OutputPixelType::EigenVectorsMatrixType EigenVectorType;
  typedef typename OutputPixelType::EigenValuesArrayType EigenValueType;
  typedef T RealType;


  /** Standard class typedefs. */
  typedef ExpEuclideanTensorImageFilter  Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  
  virtual void ThreadedGenerateData();


protected:
  ExpEuclideanTensorImageFilter() {};
  virtual ~ExpEuclideanTensorImageFilter() {};

private:
  ExpEuclideanTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkExpEuclideanTensorImageFilter.txx"
#endif

}


#endif
