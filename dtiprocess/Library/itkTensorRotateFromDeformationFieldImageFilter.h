/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorRotateFromDeformationFieldImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorRotateFromDeformationFieldImageFilter_h
#define __itkTensorRotateFromDeformationFieldImageFilter_h

#include <itkBinaryFunctorImageFilter.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_qr.h>

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput1, typename TInput2, typename TOutput >
class TensorRotateFromDeformationFieldFunction
{
public:
  typedef typename TOutput::ValueType RealType;

  TensorRotateFromDeformationFieldFunction() {}
  ~TensorRotateFromDeformationFieldFunction() {}
  bool operator!=( const TensorRotateFromDeformationFieldFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorRotateFromDeformationFieldFunction & other ) const
  {
    return !(*this != other);
  }
  TOutput operator()( const TInput1 & x, const TInput2 &y ) const
  {
    // Extract rotation from TInput2
    typedef typename TInput2::InternalMatrixType VnlMatrixType;
    typedef typename TInput2::ComponentType TransformPrecision;

    VnlMatrixType iden;
    iden.set_identity();

    VnlMatrixType vnlx;
    for(int i = 0; i < 3; ++i)
      {
      for(int j = 0; j < 3; ++j)
        {
        vnlx(i,j) = x(i,j);
        }
      }
    
    vnl_svd<TransformPrecision> svd(y.GetVnlMatrix() + iden);
    VnlMatrixType R(svd.U() * svd.V().transpose());
    
    VnlMatrixType res = R.transpose() * vnlx * R;
      
    TOutput out;
    for(int i = 0; i < 3; ++i)
      {
      for(int j = i; j < 3; ++j)
        {
        out(i,j) = res(i,j);
        }
      }

    return out;
  }

}; 


}  // end namespace functor


/** \class TensorRotateFromDeformationFieldImageFilter
 * \brief Computes the Fractional Anisotropy for every pixel of a input tensor image.
 *
 * TensorRotateFromDeformationFieldImageFilter applies pixel-wise the invokation for
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
class ITK_EXPORT TensorRotateFromDeformationFieldImageFilter :
    public
BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                        Functor::TensorRotateFromDeformationFieldFunction< 
                                        typename TInputImage1::PixelType,
                                        typename TInputImage2::PixelType,
                                        typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef TensorRotateFromDeformationFieldImageFilter  Self;
  typedef BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                                  Functor::TensorRotateFromDeformationFieldFunction< 
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
  TensorRotateFromDeformationFieldImageFilter() {};
  virtual ~TensorRotateFromDeformationFieldImageFilter() {};

private:
  TensorRotateFromDeformationFieldImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
