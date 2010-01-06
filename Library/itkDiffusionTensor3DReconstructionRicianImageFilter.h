/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionRicianImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionRicianImageFilter_h_
#define __itkDiffusionTensor3DReconstructionRicianImageFilter_h_

#include "itkImageToImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include "itkSingleValuedCostFunction.h"

#include "vnl/vnl_cost_function.h"

#include "cephes/cephes.h"

#include <iomanip>
#include <iostream>

namespace itk{

template< class TSignalType >
class RicianLikelihood : public SingleValuedCostFunction
{
public:
  typedef RicianLikelihood<TSignalType>           Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  itkNewMacro(Self);

  virtual double GetValue(const Array<double> &D) const
  {
    const double s2 = m_Sigma * m_Sigma;
    const unsigned int ns = m_Signal.size();
    
    vnl_vector<double> attenuation(ns);
    attenuation = m_Design_Matrix * D;
    for(unsigned int i = 0; i < ns; ++i)
      attenuation[i] = exp(attenuation[i]);
    vnl_vector<double> A = m_S0*attenuation;
    //int ierr;

    vnl_vector<double> loglhood(ns);
    double sumloglhood = 0.0;
    //std::cout << "s: " << m_Signal << std::endl;
    //std::cout << "A: " << A << std::endl;

//    std::cout << "s0: " << m_S0 << std::endl;
    for(unsigned int i = 0; i < ns; ++i)
      {
//      std::cout << "s2: " << s2 << std::endl;
//      std::cout << "I0: " << besseli0(m_Signal[i]*A[i]/s2,0,&ierr) << std::endl;
      //                                                                                                                                     Correction factor
//      loglhood[i] = log(m_Signal[i]) - log(s2) - (m_Signal[i]*m_Signal[i] + A[i]*A[i])/(2*s2) + log(besseli0(m_Signal[i]*A[i]/s2,1,&ierr)) + m_Signal[i]*A[i]/s2 ;
#if 0 // amos fortan bessel
      loglhood[i] =  - (m_Signal[i]*m_Signal[i] + A[i]*A[i])/(2*s2) + log(besseli0(m_Signal[i]*A[i]/s2,1,&ierr)) + m_Signal[i]*A[i]/s2 ;
#else //cephes c bessel
      loglhood[i] =  - (m_Signal[i]*m_Signal[i] + A[i]*A[i])/(2*s2) + log(i0e(m_Signal[i]*A[i]/s2)) + m_Signal[i]*A[i]/s2 ;
#endif
//       std::cout << "[" << std::setw(2) << i << "]: " << std::setw(13) << m_Signal[i]
//                 << std::setw(13) << A[i]
//                 << std::setw(13) << loglhood[i]
//                 << std::endl;
      sumloglhood += loglhood[i];
      }


//    std::cout << std::endl << "D: " << D << std::endl;
//    std::cout << "logl: " << -sumloglhood << std::endl;
    return - sumloglhood;
  }

  virtual void GetDerivative(Array<double> const& D, Array<double>& gradient) const
  {
//    std::cout << "D: "<< D << std::endl;
//    std::cout << "f(D): "<< this->GetValue(D) << std::endl;
    const double gs = 1.0e-10;
    //const double f = this->GetValue(D);
    
    Array<double> hi(6);
    Array<double> lo(6);
    for(unsigned int j = 0; j < 6; ++j)
      {
      Array<double> hid(D);
      Array<double> lod(D);
      hid[j] = hid[j] + gs;
      lod[j] = lod[j] - gs;

      hi[j] = this->GetValue(hid);
      lo[j] = this->GetValue(lod);
      gradient[j] = (hi[j] - lo[j])/(2*gs);
      }

//    std::cout << "fdg: " << gradient << std::endl;

//     const double s2 = m_Sigma * m_Sigma;
//     const unsigned int ns = m_Signal.size();
//     const vnl_matrix<double> X = m_Design_Matrix;

//     vnl_vector<double> attenuation(ns);
//     attenuation = m_Design_Matrix * D;
//     for(unsigned int i = 0; i < ns; ++i)
//       attenuation[i] = exp(attenuation[i]);
//     vnl_vector<double> A = m_S0*attenuation;
//     int ierr;

//     gradient.fill(0.0);
//     for(unsigned int j = 0; j < 6; ++j)
//       {
//         for(unsigned int i = 0; i < ns; ++i)
//         {
//         double besselarg = m_Signal[i] * m_S0 * exp(X(i,j) * D[j]) / s2;
//         gradient[j] += -(m_S0*m_S0 * X(i,j) * exp(2*X(i,j) * D[j]))/s2 + (X(i,j) * exp(X(i,j) * D[j]) * m_S0 * m_Signal[i])/s2 *
//                          besseli1(besselarg,1,&ierr) / besseli0(besselarg,1,&ierr);
//         }
//       }
//     gradient = -gradient;
//     std::cout << "g: " << gradient << std::endl;
  }

  virtual unsigned int GetNumberOfParameters() const {return 6;}

itkSetMacro(Design_Matrix, vnl_matrix<double>);
itkSetMacro(S0, double);
itkSetMacro(Sigma, double);
  
  virtual void SetSignal(const vnl_vector<TSignalType>& s)
  {
    m_Signal.set_size(s.size());
    for(unsigned i = 0; i < s.size(); ++i)
      {
      if(s[i] > 0)
        m_Signal[i] = s[i];
      else
        m_Signal[i] = 1;
      }
  }

protected:
  vnl_matrix<double> m_Design_Matrix;
  double             m_S0;
  double             m_Sigma;
  vnl_vector<double> m_Signal;

};

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType=TReferenceImagePixelType,
          class TTensorType=double >
class ITK_EXPORT DiffusionTensor3DReconstructionRicianImageFilter :
  public ImageToImageFilter< Image< TReferenceImagePixelType, 3 >, 
                             Image< DiffusionTensor3D< TTensorType >, 3 > >
{

public:

  typedef DiffusionTensor3DReconstructionRicianImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ImageToImageFilter< Image< TReferenceImagePixelType, 3>, 
          Image< DiffusionTensor3D< TTensorType >, 3 > >
                          Superclass;
  
   /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(DiffusionTensor3DReconstructionRicianImageFilter, 
                                                      ImageToImageFilter);
 
  typedef TReferenceImagePixelType                 ReferencePixelType;

  typedef TGradientImagePixelType                  GradientPixelType;

  typedef DiffusionTensor3D< TTensorType >    TensorType;

  /** Reference image data,  This image is aquired in the absence 
   * of a diffusion sensitizing field gradient */
  typedef typename Superclass::InputImageType      ReferenceImageType;
  
  typedef Image< TensorType, 3 >              TensorImageType;
  
  typedef TensorImageType                          OutputImageType;

  typedef typename Superclass::OutputImageRegionType
                                                   OutputImageRegionType;

  /** Typedef defining one (of the many) gradient images.  */
  typedef Image< GradientPixelType, 3 >            GradientImageType;

  /** An alternative typedef defining one (of the many) gradient images. 
   * It will be assumed that the vectorImage has the same dimension as the 
   * Reference image and a vector length parameter of \c n (number of
   * gradient directions)*/
  typedef VectorImage< GradientPixelType, 3 >      GradientImagesType;

  /** Holds the tensor basis coefficients G_k */
  typedef vnl_matrix< double >         TensorBasisMatrixType;
  
  typedef vnl_matrix< double >                     CoefficientMatrixType;

  /** Holds each magnetic field gradient used to acquire one DWImage */
  typedef vnl_vector_fixed< double, 3 >            GradientDirectionType;

  /** Container to hold gradient directions of the 'n' DW measurements */
  typedef VectorContainer< unsigned int, 
          GradientDirectionType >                  GradientDirectionContainerType;
  

  /** Set method to add a gradient direction and its corresponding image. */
//  void AddGradientImage( const GradientDirectionType &, const GradientImageType *image);
  
  /** Another set method to add a gradient directions and its corresponding
   * image. The image here is a VectorImage. The user is expected to pass the 
   * gradient directions in a container. The ith element of the container 
   * corresponds to the gradient direction of the ith component image the 
   * VectorImage.  For the baseline image, a vector of all zeros
   * should be set.*/
  void SetGradientImage( GradientDirectionContainerType *, 
                                              GradientImagesType *image);

  /** Set method to set initial tensor image got by linear estimation **/
  void SetInitialTensor( TensorImageType *image);

  /** Set method to set the reference image. */
  void SetReferenceImage( ReferenceImageType *referenceImage )
    {
    if( m_GradientImageTypeEnumeration == GradientIsInASingleImage)
      {
      itkExceptionMacro( << "Cannot call both methods:" 
      << "AddGradientImage and SetGradientImage. Please call only one of them.");
      }
  
    this->ProcessObject::SetNthInput( 0, referenceImage );

    m_GradientImageTypeEnumeration = GradientIsInManyImages;
    }
    
  /** Set B value to estimator **/
  //virtual void SetBValue(const TTensorType _bValue){ m_BValue = _bValue;}
   void SetStep(const TTensorType _step){ m_Step = _step;}

  /** Set Sigma **/
   void SetSigma(const TTensorType _sigma){ m_Sigma = _sigma;}
   
  /** Get reference image */
  virtual ReferenceImageType * GetReferenceImage() 
  { return ( static_cast< ReferenceImageType *>(this->ProcessObject::GetInput(0)) ); }

  /** Return the gradient direction. idx is 0 based */
  virtual GradientDirectionType GetGradientDirection( unsigned int idx) const
    {
    if( idx >= m_NumberOfGradientDirections )
      {
      itkExceptionMacro( << "Gradient direction " << idx << "does not exist" );
      }
    return m_GradientDirectionContainer->ElementAt( idx+1 );
    }

  /** Threshold on the reference image data. The output tensor will be a null
   * tensor for pixels in the reference image that have a value less than this
   * threshold. */
  itkSetMacro( Threshold, ReferencePixelType );
  itkGetMacro( Threshold, ReferencePixelType );

  
  /** 
   * The BValue \f$ (s/mm^2) \f$ value used in normalizing the tensors to 
   * physically meaningful units.  See equation (24) of the first reference for
   * a description of how this is applied to the tensor estimation.
   * Equation (1) of the same reference describes the physical significance.
   */
  itkSetMacro( BValue, TTensorType);
#ifdef GetBValue
#undef GetBValue
#endif
  itkGetConstReferenceMacro( BValue, TTensorType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(ReferenceEqualityComparableCheck,
    (Concept::EqualityComparable<ReferencePixelType>));
  itkConceptMacro(TensorEqualityComparableCheck,
    (Concept::EqualityComparable<TensorType>));
  itkConceptMacro(GradientConvertibleToDoubleCheck,
    (Concept::Convertible<GradientPixelType, double>));
  itkConceptMacro(DoubleConvertibleToTensorCheck,
    (Concept::Convertible<double, TensorType>));
  itkConceptMacro(GradientReferenceAdditiveOperatorsCheck,
    (Concept::AdditiveOperators<GradientPixelType, GradientPixelType,
                                ReferencePixelType>));
  itkConceptMacro(ReferenceOStreamWritableCheck,
    (Concept::OStreamWritable<ReferencePixelType>));
  itkConceptMacro(TensorOStreamWritableCheck,
    (Concept::OStreamWritable<TensorType>));
  /** End concept checking */
#endif

protected:
  DiffusionTensor3DReconstructionRicianImageFilter();
  ~DiffusionTensor3DReconstructionRicianImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void ComputeTensorBasis();

  //void EstimateTensor
  //double calGuessS(const TensorType & _tensor, int index);
  double ComputeGuessValue( const TensorType &_tensor, typename
         NumericTraits<ReferencePixelType>::AccumulateType cleanValue, int index );

  void BeforeThreadedGenerateData();
  void ThreadedGenerateData( const 
      OutputImageRegionType &outputRegionForThread, int);
  
  /** enum to indicate if the gradient image is specified as a single multi-
   * component image or as several separate images */
  typedef enum
    {
    GradientIsInASingleImage = 1,
    GradientIsInManyImages,
    Else
    } GradientImageTypeEnumeration;
    
private:
  
  /* Tensor basis coeffs */
  TensorBasisMatrixType                             m_TensorBasis;
  
  CoefficientMatrixType                             m_BMatrix;

  /** container to hold gradient directions */
  GradientDirectionContainerType::Pointer           m_GradientDirectionContainer;

  /** Number of gradient measurements */
  unsigned int                                      m_NumberOfGradientDirections;

  /** Number of baseline images */
  unsigned int                                      m_NumberOfBaselineImages;

  /** Threshold on the reference image data */
  ReferencePixelType                                m_Threshold;

  /** LeBihan's b-value for normalizing tensors */
  TTensorType                                  m_BValue;
 
  /** The step for SteepGradient Method **/
  TTensorType                                  m_Step;


  /** Gradient image was specified in a single image or in multiple images */
  GradientImageTypeEnumeration                      m_GradientImageTypeEnumeration;

  /** Sigma **/
  TTensorType                                  m_Sigma;

  static const double EPS = 1e-10;
  static const double ZEPS = 1e-10;


};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionRicianImageFilter.txx"
#endif

#endif

