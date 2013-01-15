/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorNegativeEigenValueImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorNegativeEigenValueImageFilter_h
#define __itkTensorNegativeEigenValueImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput, typename TOutput >
class TensorNegativeEigenValueFunction
{
public:
  TensorNegativeEigenValueFunction() {}
  ~TensorNegativeEigenValueFunction() {}
  bool operator!=( const TensorNegativeEigenValueFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorNegativeEigenValueFunction & other ) const
  {
    return !(*this != other);
  }
  TOutput operator()( const TInput & x )
    {
      
      if(x[0] == 0 && x[1] == 0 &&
         x[2] == 0 && x[3] == 0 &&
         x[4] == 0 && x[5] == 0)
        return 0;

      typedef typename TInput::EigenValuesArrayType EigenValuesType;
      EigenValuesType e;
      x.ComputeEigenValues(e);
      if(e[0] <= 0 || e[1] <= 0 || e[2] <= 0)
        return 1;
      else
        return 0;
    }
}; 

}  // end namespace functor


/** \class TensorNegativeEigenValueImageFilter
 * \brief Computes the Mean Diffusivity for every pixel of a input tensor image.
 *
 * TensorNegativeEigenValueImageFilter applies pixel-wise the invokation for
 * computing the mean diffusivity of every pixel. The pixel type of the
 * input image is expected to implement a method GetTrace(), and
 * to specify its return type as RealValueType.
 * 
 * \sa TensorRelativeAnisotropyImageFilter
 * \sa TensorFractionalAnisotropyImageFilter
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <typename TInputImage,
          typename TOutputImage>
class ITK_EXPORT TensorNegativeEigenValueImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::TensorNegativeEigenValueFunction< 
  typename TInputImage::PixelType,
  typename TOutputImage::PixelType> > 
{
public:
  /** Standard class typedefs. */
  typedef TensorNegativeEigenValueImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::TensorNegativeEigenValueFunction< 
    typename TInputImage::PixelType,
    typename TOutputImage::PixelType> >  Superclass;

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
  TensorNegativeEigenValueImageFilter() {};
  virtual ~TensorNegativeEigenValueImageFilter() {};

private:
  TensorNegativeEigenValueImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
