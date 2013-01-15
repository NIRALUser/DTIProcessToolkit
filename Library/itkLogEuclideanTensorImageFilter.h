/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLogEuclideanTensorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLogEuclideanTensorImageFilter_h
#define __itkLogEuclideanTensorImageFilter_h

#include <itkUnaryFunctorImageFilter.h>
#include <vnl/vnl_double_3x3.h>
#include <itkDiffusionTensor3D.h>

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput >
class LogEuclideanTensorFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  typedef typename itk::Vector<RealValueType, 6> OutputType;
  typedef TInput TensorType;
  typedef typename TInput::EigenVectorsMatrixType EigenVectorType;
  typedef typename TInput::EigenValuesArrayType EigenValueType;

  LogEuclideanTensorFunction() {}
  ~LogEuclideanTensorFunction() {}

  OutputType operator()( const TInput & x )
  {
    EigenValueType D;
    EigenVectorType U;

    x.ComputeEigenAnalysis(D,U);
    
    vnl_matrix_fixed<RealValueType,3,3> m;
    m.fill(0);
    m(0,0) = D[0] > 0 ? log(D[0]) : -10;
    m(1,1) = D[1] > 0 ? log(D[1]) : -10;
    m(2,2) = D[2] > 0 ? log(D[2]) : -10;

    vnl_matrix_fixed<RealValueType,3,3> matlog;
    matlog = U.GetVnlMatrix().transpose() * m * U.GetVnlMatrix();

    OutputType op;
    op[0] = matlog(0,0);
    op[1] = matlog(0,1) * sqrt(2.0);
    op[2] = matlog(0,2) * sqrt(2.0);
    op[3] = matlog(1,1);
    op[4] = matlog(1,2) * sqrt(2.0);
    op[5] = matlog(2,2);

    return op;
  }
}; 


}  // end namespace functor


/** \class LogEuclideanTensorImageFilter
 * \brief Computes the 6-element vector field that are the unique
 * elements of the matrix-logarithm of the tensor field.
 *
 * LogEuclideanImageFilter applies pixel-wise the invokation for
 * computing the matrix logarithm of every pixel. 
 * 
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template<typename T>
class ITK_EXPORT LogEuclideanTensorImageFilter :
    public
UnaryFunctorImageFilter<Image<DiffusionTensor3D<T>, 3>,
                        Image<Vector<T,6>, 3>,
                        Functor::LogEuclideanTensorFunction<DiffusionTensor3D<T> > >
{
public:
  typedef Image<DiffusionTensor3D<T>, 3> InputImageType;
  typedef Image<Vector<T, 6>, 3> OutputImageType;

  /** Standard class typedefs. */
  typedef LogEuclideanTensorImageFilter  Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
 

protected:
  LogEuclideanTensorImageFilter() {};
  virtual ~LogEuclideanTensorImageFilter() {};

private:
  LogEuclideanTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



}


#endif
