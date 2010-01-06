/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorFrobeniusNormImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorFrobeniusNormImageFilter_h
#define __itkTensorFrobeniusNormImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{

// This functor class invokes the computation of mean diffusivity from
// every pixel.
namespace Functor {  
 
template< typename TInput >
class TensorFrobeniusNormFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  TensorFrobeniusNormFunction() {}
  ~TensorFrobeniusNormFunction() {}
  bool operator!=( const TensorFrobeniusNormFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorFrobeniusNormFunction & other ) const
  {
    return !(*this != other);
  }
  inline RealValueType operator()( const TInput & x )
    {
      return sqrt(x[0]*x[0] +
                  2*x[1]*x[1] +
                  2*x[2]*x[2] +
                  x[3]*x[3] + 
                  2*x[4]*x[4] + 
                  x[5]*x[5]);
    }
}; 

}  // end namespace functor


/** \class TensorFrobeniusNormImageFilter
 * \brief Computes the Mean Diffusivity for every pixel of a input tensor image.
 *
 * TensorFrobeniusNormImageFilter applies pixel-wise the invokation for
 * computing the mean diffusivity of every pixel. The pixel type of the
 * input image is expected to implement a method GetTrace(), and
 * to specify its return type as RealValueType.  The mean diffusivity
 * is the average diffusion over all directions.
 * 
 * \sa TensorRelativeAnisotropyImageFilter
 * \sa TensorFractionalAnisotropyImageFilter
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <typename TInputImage,
          typename TOutputImage=itk::Image<ITK_TYPENAME TInputImage::PixelType::RealValueType,
                                           ::itk::GetImageDimension<TInputImage>::ImageDimension > >
class ITK_EXPORT TensorFrobeniusNormImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::TensorFrobeniusNormFunction< 
                                        typename TInputImage::PixelType> > 
{
public:
  /** Standard class typedefs. */
  typedef TensorFrobeniusNormImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::TensorFrobeniusNormFunction< 
                                    typename TInputImage::PixelType> >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename TInputImage::PixelType         InputPixelType;
  typedef typename InputPixelType::ValueType      InputValueType;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  

protected:
  TensorFrobeniusNormImageFilter() {};
  virtual ~TensorFrobeniusNormImageFilter() {};

private:
  TensorFrobeniusNormImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
