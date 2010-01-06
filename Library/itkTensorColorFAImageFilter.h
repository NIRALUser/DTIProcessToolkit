/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorColorFAImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorColorFAImageFilter_h
#define __itkTensorColorFAImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkRGBPixel.h"

namespace itk
{

// This functor class invokes the eigendecomposition of a tensor
// pixel.  The returned value is the major eigenvector in RGB space
// weighted by the FA value.
namespace Functor {  
 
template< typename TInput , typename RGBPixelComponentType>
class TensorColorFAFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  typedef typename TInput::EigenVectorsMatrixType EigenVectorsType;
  typedef typename TInput::EigenValuesArrayType EigenValuesType;
  typedef RGBPixel<RGBPixelComponentType> PixelType;

  TensorColorFAFunction() {}
  ~TensorColorFAFunction() {}
  bool operator!=( const TensorColorFAFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorColorFAFunction & other ) const
  {
    return !(*this != other);
  }
  PixelType operator()( const TInput & x )
    {
      RealValueType fa = x.GetFractionalAnisotropy();
      // Clamp FA
      if(fa > 1.0) fa = 1.0;

      PixelType color;
      EigenVectorsType mat;
      EigenValuesType e;
      RealValueType ev1[3];

      x.ComputeEigenAnalysis(e,mat);

      if(e[1] > e[0] && e[1] > e[2])
        {
        ev1[0] = mat(1,0); ev1[1] = mat(1,1); ev1[2] = mat(1,2);
        }
      else if(e[2] > e[0] && e[2] > e[1])
        {
        ev1[0] = mat(2,0); ev1[1] = mat(2,1); ev1[2] = mat(2,2);
        }
      else
        {
        ev1[0] = mat(0,0); ev1[1] = mat(0,1); ev1[2] = mat(0,2);
        }
      
      color.Set(static_cast<RGBPixelComponentType>(fabs(ev1[0]) * fa * NumericTraits<RGBPixelComponentType>::max()),
                static_cast<RGBPixelComponentType>(fabs(ev1[1]) * fa * NumericTraits<RGBPixelComponentType>::max()),
                static_cast<RGBPixelComponentType>(fabs(ev1[2]) * fa * NumericTraits<RGBPixelComponentType>::max()));
      return color;
    }
}; 

}  // end namespace functor


/** \class TensorColorFAImageFilter
 * \brief Computes the color FA value for every pixel of an input tensor image.
 *
 * TensorColorFAImageFilter applies pixel-wise the invokation for
 * computing the eigenvectors and fractional anisotropy of a tensor
 * pixel. The pixel type of the input image is expected to implement a
 * method GetFractionalAnisotropy() and ComputeEigenAnalysis().  The
 * resulting color is the eigenvector times the FA value.
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
class ITK_EXPORT TensorColorFAImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::TensorColorFAFunction<
                              typename TInputImage::PixelType,
                              typename TOutputImage::PixelType::ComponentType> > 
{
public:
  /** Standard class typedefs. */
  typedef TensorColorFAImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::TensorColorFAFunction< 
                                    typename TInputImage::PixelType,
                                    typename TOutputImage::PixelType::ComponentType> >  Superclass;

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
  TensorColorFAImageFilter() {};
  virtual ~TensorColorFAImageFilter() {};

private:
  TensorColorFAImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
