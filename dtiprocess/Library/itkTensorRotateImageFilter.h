/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorRotateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorRotateImageFilter_h
#define __itkTensorRotateImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput, typename TOutput, typename TTransformRealType >
class TensorRotateFunction
{
public:
  typedef TTransformRealType RealType;

  TensorRotateFunction() {}
  ~TensorRotateFunction() {}
  bool operator!=( const TensorRotateFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorRotateFunction & other ) const
  {
    return !(*this != other);
  }
  TOutput operator()( const TInput & x )
  {
    vnl_matrix<RealType> T(3,3);
    for(int i = 0; i < 3; ++i)
      {
      for(int j = 0; j < 3; ++j)
        {
        T(i,j) = x(i,j);
        }
      }

    // TODO is the transpose in the right order?
    vnl_matrix<RealType> Tr = m_Rotation * T * m_Rotation.transpose();

    TOutput ret;
    for(int i = 0; i < 3; ++i)
      {
      for(int j = 0; j < 3; ++j)
        {
        ret(i,j) = Tr(i,j);
        }
      }
    return ret;
  }

  vnl_matrix<RealType> GetRotation() const
  {
    return m_Rotation;
  }

  void SetRotation(const vnl_matrix<RealType> &a)
  {
    m_Rotation = a;
  }

private:
  vnl_matrix<RealType> m_Rotation;
}; 

}  // end namespace functor


/** \class TensorRotateImageFilter
 * \brief Computes the Fractional Anisotropy for every pixel of a input tensor image.
 *
 * TensorRotateImageFilter applies pixel-wise the invokation for
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
template <typename TInputImage,
          typename TOutputImage,
          typename TTransformRealType>
class ITK_EXPORT TensorRotateImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::TensorRotateFunction< 
                                        typename TInputImage::PixelType,
                                        typename TOutputImage::PixelType,
                                        TTransformRealType> > 
{
public:
  /** Standard class typedefs. */
  typedef TensorRotateImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::TensorRotateFunction< 
                                    typename TInputImage::PixelType,
                                    typename TOutputImage::PixelType,
                                    TTransformRealType> >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename TInputImage::PixelType         InputPixelType;
  typedef typename InputPixelType::ValueType      InputValueType;

  typedef TTransformRealType TransformRealType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }

  vnl_matrix<TransformRealType> GetRotation() const
  {
    return this->GetFunctor().GetRotation();
  }

  void SetRotation(const vnl_matrix<TransformRealType> &a)
  {
    this->GetFunctor().SetRotation(a);
    this->Modified();
  }

protected:
  TensorRotateImageFilter() {};
  virtual ~TensorRotateImageFilter() {};

private:
  TensorRotateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
