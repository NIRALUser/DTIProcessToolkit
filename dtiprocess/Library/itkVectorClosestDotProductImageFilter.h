/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorClosestDotProductImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorClosestDotProductImageFilter_h
#define __itkVectorClosestDotProductImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkCovariantVector.h"
#include "itkVectorContainer.h"

namespace itk
{

// This functor class invokes the computation of fractional anisotropy from
// every pixel.
namespace Functor {  
 
template< typename TInput , typename TOutput>
class VectorClosestDotProductFunction
{
public:
  typedef vnl_vector_fixed<double, 3> GradientType;
  typedef VectorContainer<unsigned int, GradientType> GradientListType;

  VectorClosestDotProductFunction() {}
  ~VectorClosestDotProductFunction() {}
  bool operator!=( const VectorClosestDotProductFunction & ) const
  {
    return false;
  }
  bool operator==( const VectorClosestDotProductFunction & other ) const
  {
    return !(*this != other);
  }
  TOutput operator()( const TInput & x )
  {
    TOutput max = NumericTraits<TOutput>::min();

    for(unsigned int i = 0; i < m_GradientList->Size(); ++i)
      {
      TOutput dot = fabs(dot_product(x.GetVnlVector(), m_GradientList->ElementAt(i)));
      if(dot > max)
        max = dot;
      }
    
    return max;
  }

  void SetGradientList(typename GradientListType::Pointer g)
  {
    m_GradientList = GradientListType::New();

    unsigned int ind = 0;
    for(unsigned int i = 0; i < g->Size(); ++i)
      {
      if(g->ElementAt(i).one_norm() != 0.0)
        {
        m_GradientList->InsertElement(ind,g->ElementAt(i));
        ind++;
        }
      }
  }

private:
  typename GradientListType::Pointer m_GradientList;

}; 

}  // end namespace functor


/** \class VectorClosestDotProductImageFilter
 * \brief Computes the Mean Diffusivity for every pixel of a input tensor image.
 *
 * VectorClosestDotProductImageFilter applies pixel-wise the invokation for
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
class ITK_EXPORT VectorClosestDotProductImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::VectorClosestDotProductFunction<
                              typename TInputImage::PixelType,
                              typename TOutputImage::PixelType> > 
{
public:
  /** Standard class typedefs. */
  typedef VectorClosestDotProductImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::VectorClosestDotProductFunction< 
                                    typename TInputImage::PixelType,
                                    typename TOutputImage::PixelType> >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename TInputImage::PixelType         InputPixelType;

  typedef typename Functor::VectorClosestDotProductFunction< 
    typename TInputImage::PixelType,
    typename TOutputImage::PixelType>  FunctorType;

  typedef typename FunctorType::GradientType GradientType;
  typedef typename FunctorType::GradientListType GradientListType;
  typedef typename GradientListType::Pointer GradientListPointerType;

  void SetGradientList(GradientListPointerType g)
  {
    this->GetFunctor().SetGradientList(g);
    this->Modified();
  }

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  

protected:
  VectorClosestDotProductImageFilter() {};
  virtual ~VectorClosestDotProductImageFilter() {};

private:
  VectorClosestDotProductImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
