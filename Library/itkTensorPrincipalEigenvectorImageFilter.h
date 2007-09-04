/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorPrincipalEigenvectorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007-09-04 20:12:29 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorPrincipalEigenvectorImageFilter_h
#define __itkTensorPrincipalEigenvectorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkCovariantVector.h"

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput , typename VectorPixelValueType>
class TensorPrincipalEigenvectorFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  typedef typename TInput::EigenVectorsMatrixType EigenVectorsType;
  typedef typename TInput::EigenValuesArrayType EigenValuesType;
  typedef CovariantVector<VectorPixelValueType,3> PixelType;

  TensorPrincipalEigenvectorFunction() {}
  ~TensorPrincipalEigenvectorFunction() {}
  bool operator!=( const TensorPrincipalEigenvectorFunction & ) const
  {
    return false;
  }
  bool operator==( const TensorPrincipalEigenvectorFunction & other ) const
  {
    return !(*this != other);
  }
  PixelType operator()( const TInput & x )
    {
      EigenVectorsType mat;
      EigenValuesType e;

      x.ComputeEigenAnalysis(e,mat);

      PixelType ev1;
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

      return ev1;
    }
}; 

}  // end namespace functor


/** \class TensorPrincipalEigenvectorImageFilter
 * \brief Computes the Mean Diffusivity for every pixel of a input tensor image.
 *
 * TensorPrincipalEigenvectorImageFilter applies pixel-wise the invokation for
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
class ITK_EXPORT TensorPrincipalEigenvectorImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::TensorPrincipalEigenvectorFunction<
                              typename TInputImage::PixelType,
                              typename TOutputImage::PixelType::ValueType> > 
{
public:
  /** Standard class typedefs. */
  typedef TensorPrincipalEigenvectorImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::TensorPrincipalEigenvectorFunction< 
                                    typename TInputImage::PixelType,
                                    typename TOutputImage::PixelType::ValueType> >  Superclass;

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
  TensorPrincipalEigenvectorImageFilter() {};
  virtual ~TensorPrincipalEigenvectorImageFilter() {};

private:
  TensorPrincipalEigenvectorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
