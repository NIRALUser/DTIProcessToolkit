/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorRotateFromDeformationFieldPPDImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorRotateFromDeformationFieldPPDImageFilter_h
#define __itkTensorRotateFromDeformationFieldPPDImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include <vnl/algo/vnl_qr.h>
#include <vnl/vnl_quaternion.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_inverse.h>

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput1, typename TInput2, typename TOutput >
class TensorRotateFromDeformationFieldPPDFunction
{
public:
  typedef typename TOutput::ValueType RealType;

  TensorRotateFromDeformationFieldPPDFunction() {}
  ~TensorRotateFromDeformationFieldPPDFunction() {}
  bool operator!=( const TensorRotateFromDeformationFieldPPDFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorRotateFromDeformationFieldPPDFunction & other ) const
  {
    return !(*this != other);
  }
  TOutput operator()( const TInput1 & x, const TInput2 &y )
  {
    // Extract rotation from TInput2
    typedef typename TInput2::InternalMatrixType VnlMatrixType;
    typedef typename TInput2::ComponentType TransformPrecision;
    VnlMatrixType iden;
    iden.set_identity();

    // Add jacobian to identity and invert.  This requires that 
    // the deformation field be invertible locally at every point
    // otherwise A will not be invertible.
    VnlMatrixType A(y.GetVnlMatrix() + iden);
    A = vnl_inverse(A);

    typedef typename TInput1::RealValueType  RealValueType;
    typedef typename TInput1::EigenVectorsMatrixType EigenVectorsType;
    typedef typename TInput1::EigenValuesArrayType EigenValuesType;
    typedef CovariantVector<RealValueType,3> PixelType;

    EigenVectorsType mat;
    EigenValuesType e;
    
    x.ComputeEigenAnalysis(e,mat);

    typedef vnl_vector_fixed<TransformPrecision, 3> VnlVectorType;
    VnlVectorType ev1, ev2, ev3;
    VnlVectorType n1, n2, pn2;
    typedef vnl_quaternion<TransformPrecision> RotationType;

    // find largest eigenvector
    ev1[0] = mat(2,0); ev1[1] = mat(2,1); ev1[2] = mat(2,2);
    ev2[0] = mat(1,0); ev2[1] = mat(1,1); ev2[2] = mat(1,2);
    ev3[0] = mat(0,0); ev3[1] = mat(0,1); ev3[2] = mat(0,2);

    n1 = (A * ev1).normalize();
    RotationType R1(vnl_cross_3d(ev1,n1).normalize(),angle(n1,ev1));

    n2 = (A * ev2).normalize();
    pn2 = n2 - dot_product(n1,n2)*n1;
    RotationType R2(R1.rotate(ev1),dot_product(R1.rotate(ev2),pn2.normalize()));

    VnlMatrixType R = (R2 * R1).rotation_matrix_transpose();

    VnlMatrixType vnlx;
    for(int i = 0; i < 3; ++i)
      {
      for(int j = 0; j < 3; ++j)
        {
        vnlx(i,j) = x(i,j);
        }
      }

     VnlMatrixType res = R * vnlx * R.transpose();
     
     TOutput out;
     for(int i = 0; i < 3; ++i)
     {
       for(int j = 0; j < 3; ++j)
       {
         out(i,j) = res(i,j);
       }
     }
     
     return out;
  }
}; 

}  // end namespace functor


/** \class TensorRotateFromDeformationFieldPPDImageFilter
 * \brief Computes the Fractional Anisotropy for every pixel of a input tensor image.
 *
 * TensorRotateFromDeformationFieldPPDImageFilter applies pixel-wise the invokation for
 * computing the fractional anisotropy of every pixel. The pixel type of the
 * input image is expected to implement a method GetFractionalAnisotropy(), and
 * to specify its return type as  RealValueType.
 * 
 * \sa TensorRelativeAnisotropyImageFilter
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <typename TInputImage1,
          typename TInputImage2,
          typename TOutputImage>
class ITK_EXPORT TensorRotateFromDeformationFieldPPDImageFilter :
    public
BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                        Functor::TensorRotateFromDeformationFieldPPDFunction< 
                                        typename TInputImage1::PixelType,
                                        typename TInputImage2::PixelType,
                                        typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef TensorRotateFromDeformationFieldPPDImageFilter  Self;
  typedef BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                                  Functor::TensorRotateFromDeformationFieldPPDFunction< 
                                    typename TInputImage1::PixelType,
                                    typename TInputImage2::PixelType,
                                    typename TOutputImage::PixelType> >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename TInputImage1::PixelType         InputPixelType1;
  typedef typename InputPixelType1::ValueType      InputValueType1;
  typedef typename TInputImage2::PixelType         InputPixelType2;
  typedef typename InputPixelType2::ValueType      InputValueType2;

  typedef typename OutputPixelType::ValueType TransformRealType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }

protected:
  TensorRotateFromDeformationFieldPPDImageFilter() {};
  virtual ~TensorRotateFromDeformationFieldPPDImageFilter() {};

private:
  TensorRotateFromDeformationFieldPPDImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
