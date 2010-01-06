/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFastSymmetricEigenAnalysisImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSymmetricEigenAnalysisImageFilter_h
#define __itkSymmetricEigenAnalysisImageFilter_h

#include <itkUnaryFunctorImageFilter.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace itk
{

// This functor class invokes the computation of Eigen Analysis for
// every pixel. The input pixel type must provide the API for the [][]
// operator, while the output pixel type must provide the API for the
// [] operator. Input pixel matrices should be symmetric.
// 
// The default operation is to order eigen values in ascending order.
// You may also use OrderEigenValuesBy( ) to order eigen values by
// magnitude as is common with use of tensors in vessel extraction.
namespace Functor {  
 
template< typename TInput, typename TOutput >
class SymmetricEigenAnalysisFunction
{
public:
  typedef typename TInput::RealValueType  RealValueType;
  SymmetricEigenAnalysisFunction() {}
  ~SymmetricEigenAnalysisFunction() {}

  bool operator!=( const SymmetricEigenAnalysisFunction & ) const
    {
    return false;
    }
  bool operator==( const SymmetricEigenAnalysisFunction & other ) const
    {
    return !(*this != other);
    }

  inline TOutput operator()( const TInput & x )
    {
    double lambdas[3];

    vnl_symmetric_eigensystem_compute_eigenvals(x[0], x[1], x[2],
                                                      x[3], x[4],
                                                            x[5],
                                                lambdas[0],
                                                lambdas[1],
                                                lambdas[2]);   
    return TOutput(lambdas);
    }

  /** Typdedefs to order eigen values. 
   * OrderByValue:      lambda_1 < lambda_2 < ....
   * OrderByMagnitude:  |lambda_1| < |lambda_2| < .....
   * DoNotOrder:        Default order of eigen values obtained after QL method
   */
  typedef enum {
    OrderByValue=1,
    OrderByMagnitude,
    DoNotOrder
  } EigenValueOrderType;
 
  /** Order eigen values. Default is to OrderByValue:  lambda_1 <
   * lambda_2 < .... */
  void OrderEigenValuesBy( EigenValueOrderType order )
    {
      m_Order = order;
    }

private:
  EigenValueOrderType m_Order;
}; 

}  // end namespace functor


/** \class SymmetricEigenAnalysisImageFilter
 * \brief Computes the Fractional Anisotropy for every pixel of a input tensor image.
 *
 * SymmetricEigenAnalysisImageFilter applies pixel-wise the invokation for
 * computing the fractional anisotropy of every pixel. The pixel type of the
 * input image is expected to implement a method GetFractionalAnisotropy(), and
 * to specify its return type as  RealValueType.
 * 
 * The OrderEigenValuesBy( .. ) method can be used to order eigen values 
 * in ascending order by value or magnitude or no ordering.
 * OrderByValue:      lambda_1 < lambda_2 < ....
 * OrderByMagnitude:  |lambda_1| < |lambda_2| < .....
 * DoNotOrder:        Default order of eigen values obtained after QL method
 *
 * The user of this class is explicitly supposed to set the dimension of the 
 * 2D matrix using the SetDimension() method.
 *
 * \sa TensorRelativeAnisotropyImageFilter
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <typename  TInputImage, typename  TOutputImage=TInputImage>
class ITK_EXPORT FastSymmetricEigenAnalysisImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::SymmetricEigenAnalysisFunction< 
                                        typename TInputImage::PixelType,
                                        typename TOutputImage::PixelType > >
{
public:
  /** Standard class typedefs. */
  typedef FastSymmetricEigenAnalysisImageFilter  Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Functor::SymmetricEigenAnalysisFunction< 
      typename TInputImage::PixelType,
      typename TOutputImage::PixelType > >   Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef SmartPointer<const Self>           ConstPointer;

  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef typename TOutputImage::PixelType        OutputPixelType;
  typedef typename TInputImage::PixelType         InputPixelType;
  typedef typename InputPixelType::ValueType      InputValueType;
  typedef typename Superclass::FunctorType        FunctorType; 

  /** Typdedefs to order eigen values. 
   * OrderByValue:      lambda_1 < lambda_2 < ....
   * OrderByMagnitude:  |lambda_1| < |lambda_2| < .....
   * DoNotOrder:        Default order of eigen values obtained after QL method
   */
  typedef typename FunctorType::EigenValueOrderType         EigenValueOrderType;
 
  /** Order eigen values. Default is to OrderByValue:  lambda_1 <
   * lambda_2 < .... */
  void OrderEigenValuesBy( EigenValueOrderType order )
    {
    this->GetFunctor().OrderEigenValuesBy( order );
    }

  /** Run-time type information (and related methods).   */
  itkTypeMacro( FastSymmetricEigenAnalysisImageFilter, UnaryFunctorImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<InputValueType>));
  /** End concept checking */
#endif

protected:
  FastSymmetricEigenAnalysisImageFilter() {};
  virtual ~FastSymmetricEigenAnalysisImageFilter() {};

private:
  FastSymmetricEigenAnalysisImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


  
} // end namespace itk
  
#endif
