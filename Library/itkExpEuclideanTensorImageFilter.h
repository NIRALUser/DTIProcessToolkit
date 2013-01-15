/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkExpEuclideanTensorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkExpEuclideanTensorImageFilter_h
#define __itkExpEuclideanTensorImageFilter_h

#include <itkUnaryFunctorImageFilter.h>
#include <vnl/vnl_double_3x3.h>
#include <itkDiffusionTensor3D.h>

namespace itk
{

//This functor class invokes the computation of fractional anisotropy from
//every pixel.
namespace Functor {  
 
template< typename TInput >
class ExpEuclideanTensorFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  typedef typename itk::DiffusionTensor3D<RealValueType> OutputType;
  typedef OutputType TensorType;
  typedef typename OutputType::EigenVectorsMatrixType EigenVectorType;
  typedef typename OutputType::EigenValuesArrayType EigenValueType;

  ExpEuclideanTensorFunction() {}
  ~ExpEuclideanTensorFunction() {}
  OutputType operator()( const TInput & x )
  {
    DiffusionTensor3D<double> tensor;
    for(int i = 0; i < 6; ++i)
      tensor[i] = x[i];
    tensor[1] /= sqrt(2.0);
    tensor[2] /= sqrt(2.0);
    tensor[4] /= sqrt(2.0);
    
    EigenValueType D;
    EigenVectorType U;

    tensor.ComputeEigenAnalysis(D,U);

    vnl_matrix_fixed<double,3,3> m;
    m.fill(0);
    m(0,0) = exp(D[0]);
    m(1,1) = exp(D[1]);
    m(2,2) = exp(D[2]);
    
//    std::cout << U << std::endl;
//    std::cout << D << std::endl;

    vnl_matrix_fixed<double,3,3> res(U.GetVnlMatrix().transpose() * m * U.GetVnlMatrix());

    OutputType op;
    op[0] = res(0,0);
    op[1] = res(0,1);
    op[2] = res(0,2);
    op[3] = res(1,1);
    op[4] = res(1,2);
    op[5] = res(2,2);

    return op;
  }
}; 


}  // end namespace functor


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
UnaryFunctorImageFilter<Image<Vector<T,6>, 3>,
                   Image<DiffusionTensor3D<T>, 3>,
                   Functor::ExpEuclideanTensorFunction<Vector<T, 6> > >
{
public:
//  typedef Image<Vector<typename TInputImage::PixelType::RealValueType,6>,3 > TOutputImage;
  typedef Image<Vector<T,6>, 3> InputImageType;
  typedef Image<DiffusionTensor3D<T>, 3> OutputImageType;

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
  
protected:
  ExpEuclideanTensorImageFilter() {};
  virtual ~ExpEuclideanTensorImageFilter() {};

private:
  ExpEuclideanTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



}


#endif
