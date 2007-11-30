/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionWeightedImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007-11-30 18:44:14 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionWeightedImageFilter_h_
#define __itkDiffusionTensor3DReconstructionWeightedImageFilter_h_

#include "itkImageToImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_copy.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"

namespace itk{

template< class TReferenceImagePixelType, 
          class TGradientImagePixelType=TReferenceImagePixelType,
          class TTensorType=double >
class ITK_EXPORT DiffusionTensor3DReconstructionWeightedImageFilter :
  public ImageToImageFilter< Image< TReferenceImagePixelType, 3 >, 
                             Image< DiffusionTensor3D< TTensorType >, 3 > >
{

public:

  typedef DiffusionTensor3DReconstructionWeightedImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ImageToImageFilter< Image< TReferenceImagePixelType, 3>, 
          Image< DiffusionTensor3D< TTensorType >, 3 > >
                          Superclass;
  
   /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(DiffusionTensor3DReconstructionWeightedImageFilter, 
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
  typedef vnl_matrix< double >                     TensorBasisMatrixType;
  
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

  itkSetMacro( NumberOfIterations, unsigned int );
  itkGetMacro( NumberOfIterations, unsigned int );
  
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
  DiffusionTensor3DReconstructionWeightedImageFilter();
  ~DiffusionTensor3DReconstructionWeightedImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void ComputeTensorBasis();
  //void EstimateTensor

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
  GradientImageTypeEnumeration                 m_GradientImageTypeEnumeration;
  
  /** */
  unsigned int                                 m_NumberOfIterations;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.txx"
#endif

#endif

