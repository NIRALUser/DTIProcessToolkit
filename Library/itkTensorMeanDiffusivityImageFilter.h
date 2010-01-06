/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorMeanDiffusivityImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorMeanDiffusivityImageFilter_h
#define __itkTensorMeanDiffusivityImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{

// This functor class invokes the computation of mean diffusivity from
// every pixel.
namespace Functor {  
 
template< typename TInput >
class TensorMeanDiffusivityFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  TensorMeanDiffusivityFunction() {}
  ~TensorMeanDiffusivityFunction() {}
  bool operator!=( const TensorMeanDiffusivityFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorMeanDiffusivityFunction & other ) const
  {
    return !(*this != other);
  }
  inline RealValueType operator()( const TInput & x )
    {
      return x.GetTrace() / 3.0;
    }
}; 

}  // end namespace functor


/** \class TensorMeanDiffusivityImageFilter
 * \brief Computes the Mean Diffusivity for every pixel of a input tensor image.
 *
 * TensorMeanDiffusivityImageFilter applies pixel-wise the invokation for
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
class ITK_EXPORT TensorMeanDiffusivityImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::TensorMeanDiffusivityFunction< 
                                        typename TInputImage::PixelType> > 
{
public:
  /** Standard class typedefs. */
  typedef TensorMeanDiffusivityImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::TensorMeanDiffusivityFunction< 
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
  TensorMeanDiffusivityImageFilter() {};
  virtual ~TensorMeanDiffusivityImageFilter() {};

private:
  TensorMeanDiffusivityImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
