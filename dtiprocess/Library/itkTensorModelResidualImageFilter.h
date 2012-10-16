/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorModelResidualImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorModelResidualImageFilter_h
#define __itkTensorModelResidualImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"


namespace itk
{
  
/** \class TensorModelResidualImageFilter
 * \brief Implements an operator for pixel-wise masking of the input 
 * image with the negative of a mask.
 *
 * This class is parametrized over the types of the  
 * input image type, the mask image type and the type of the output image. 
 * Numeric conversions (castings) are done by the C++ defaults.
 *
 * The pixel type of the input 2 image must have a valid defintion of the
 * operator != with zero . This condition is required because internally this
 * filter will perform the operation
 *
 *        if pixel_from_mask_image != 0 
 *             pixel_output_image = 0
 *        else
 *             pixel_output_image = pixel_input_image
 *
 * The pixel from the input 1 is cast to the pixel type of the output image.
 *
 * Note that the input and the mask images must be of the same size.
 *
 * \warning Any pixel value other than 0 will not be masked out. 
 *
 * \sa MaskImageFilter
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  
  
template< class TDWIVector, class TTensor, class TTensorOut>
class TensorModelResidualInput
{
public:
  typedef vnl_vector_fixed<double, 3> GradientType;
  typedef VectorContainer<unsigned int, GradientType> GradientListType;
  typedef typename TDWIVector::ComponentType DWIComponentType;
  typedef typename NumericTraits< DWIComponentType >::AccumulateType AccumulatorType;

  TensorModelResidualInput() {};
  ~TensorModelResidualInput() {};

  bool operator!=( const TensorModelResidualInput & ) const
  {
    return false;
  }

  bool operator==( const TensorModelResidualInput & other ) const
  {
    return !(*this != other);
  }

  TTensor operator()( const TDWIVector & A, const TTensor & D)
  {
    AccumulatorType baseline = NumericTraits<AccumulatorType>::Zero;
    VariableLengthVector<double> proj(A.Size() - m_NumBaselines);
    VariableLengthVector<double> trunc(A.Size() - m_NumBaselines);

    TTensor Dout;

    unsigned int gradnum = 0;
    for(unsigned int i = 0; i < A.Size(); ++i)
      {
      if(m_GradientList->ElementAt(i).one_norm() == 0.0)
        {
        baseline += A[i];
        }
      else
        {
        trunc[gradnum] = A[i];
        gradnum++;
        }
      }
    
    assert(gradnum + m_NumBaselines == A.Size());

    baseline /= static_cast<AccumulatorType>(m_NumBaselines);

    TResidual residual = NumericTraits<AccumulatorType>::Zero;

    // ignore this pixel if a tensor was not estimated
    if(D[0] == 0 &&
       D[1] == 0 &&
       D[2] == 0 &&
       D[3] == 0 &&
       D[4] == 0 &&
       D[5] == 0)

      return 0;

    gradnum = 0;
    for(unsigned int i = 0; i < A.Size(); ++i)
      {
      if(m_GradientList->ElementAt(i).one_norm() != 0.0)
        {
        GradientType g = m_GradientList->ElementAt(i);

        double gdg = g[0]*(D(0,0)*g[0] + D(0,1)*g[1] + D(0,2)*g[2])+
                     g[1]*(D(1,0)*g[0] + D(1,1)*g[1] + D(1,2)*g[2])+
                     g[2]*(D(2,0)*g[0] + D(2,1)*g[1] + D(2,2)*g[2]);

        proj[gradnum] = baseline * exp(-m_BValue * gdg);

        gradnum++;
        }
      
      }
    
    return Dout;
  }

  void SetGradientList(typename GradientListType::Pointer g)
  {
    m_GradientList = g;

    m_NumBaselines = 0;
    for(unsigned int i = 0; i < m_GradientList->Size(); ++i)
      {
      if(m_GradientList->ElementAt(i).one_norm() == 0.0)
        {
        m_NumBaselines++;
        }
      }
    assert(m_NumBaselines > 0);
  }

  void SetBValue(double b)
  {
    m_BValue = b;
  }

private:
  double m_BValue;
  typename GradientListType::Pointer m_GradientList;
  unsigned int m_NumBaselines;
  
}; 

}

template <class TDWIImage, class TTensorImage, class TTensorOutImage>
class ITK_EXPORT TensorModelResidualImageFilter :
    public
    BinaryFunctorImageFilter<TDWIImage,TTensorImage,TTensorOutImage, 
                         Functor::TensorModelResidualInput< 
      typename TDWIImage::PixelType, 
      typename TTensorImage::PixelType,
      typename TTensorOutImage::PixelType>   >


{
public:
  /** Standard class typedefs. */
  typedef TensorModelResidualImageFilter  Self;
  typedef typename Functor::TensorModelResidualInput< 
    typename TDWIImage::PixelType, 
    typename TTensorImage::PixelType,
    typename TTensorOutImage::PixelType>   FunctorType;

  typedef BinaryFunctorImageFilter<TDWIImage,
                                    TTensorImage,
                                    TTensorOutImage, 
                                    FunctorType>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Typedefs for Functor */
  typedef typename FunctorType::GradientType GradientType;
  typedef typename FunctorType::GradientListType GradientListType;
  typedef typename GradientListType::Pointer GradientListPointerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  void SetBValue(double b)
  {
    this->GetFunctor().SetBValue(b);
    this->Modified();
  }

  void SetGradientList(GradientListPointerType g)
  {
    this->GetFunctor().SetGradientList(g);
    this->Modified();
  }

protected:
  TensorModelResidualImageFilter() {}
  virtual ~TensorModelResidualImageFilter() {}

private:
  TensorModelResidualImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
