/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.3 $

  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionImageFilterBase_h_
#define __itkDiffusionTensor3DReconstructionImageFilterBase_h_

#include "itkImageToImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"

namespace itk{
/** \class DiffusionTensor3DReconstructionImageFilterBase 
 *
 *\brief This class is a base class for filters that take as input
 * vector images with each vector representing a set of diffusion
 * weighted signals, their gradient directions and computes an image
 * of tensors. (with DiffusionTensor3D as the pixel type). 
 * 
 * Once that is done, you can apply filters on this tensor image to
 * compute FA, ADC, RGB weighted maps etc.
 *
 * \par Inputs and Usage
 *
 * The diffusion weighted images must be provided as a vector image
 * with a container listing the vector images.
 * \code
 *       filter->SetGradientImage( directionsContainer, vectorImage );
 * \endcode
 * This is convenient when the DWI images are read in using the 
 * <a href="http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format">NRRD</a> 
 * format. Like the Nrrd format, the reference images are those components of the 
 * vectorImage whose gradient direction is (0,0,0).  The magnitude of
 * the vectors scales the b-value by the magnitude of gradient
 * vector.  This a vector (2,0,0) with b-value 1000 is the same as the
 * vector (1,0,0) with b-value of 2000.
 *
 * \par Outputs
 * The output image is an image of Tensors:
 * \code
 *       Image< DiffusionTensor3D< TTensorPrecision >, 3 >
 * \endcode
 *
 * \par Parameters
 * \li Threshold -  Threshold on the reference image data. The output tensor will 
 * be a null tensor for pixels in the reference image that have a value less 
 * than this.
 * \li BValue - See the documentation of SetBValue().  This is the
 * b-value for gradients with unit norm.
 * \li At least 7 gradient images must be specified for the filter to be able 
 * to run.
 * 
 * \par Template parameters
 * The class is templated over the pixel type of the reference and gradient 
 * images (expected to be scalar data types) and the internal representation
 * of the DiffusionTensor3D pixel (double, float etc).
 *  
 * \par References:
 * \li<a href="http://lmi.bwh.harvard.edu/papers/pdfs/2002/westinMEDIA02.pdf">[1]</a> 
 * <em>C.F.Westin, S.E.Maier, H.Mamata, A.Nabavi, F.A.Jolesz, R.Kikinis,
 * "Processing and visualization for Diffusion tensor MRI", Medical image
 * Analysis, 2002, pp 93-108.</em>
 * \li<a href="splweb.bwh.harvard.edu:8000/pages/papers/westin/ISMRM2002.pdf">[2]</a>
 * <em>A Dual Tensor Basis Solution to the Stejskal-Tanner Equations for DT-MRI</em>
 * 
 * \author Thanks to Xiaodong Tao, GE, for contributing the original
 * version of this class.
 * 
 * \note
 * This work is part of the National Alliance for Medical image Computing 
 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
 * for Medical Research, Grant U54 EB005149.
 *
 * \par Examples and Datasets
 * See Examples/Filtering/DiffusionTensor3DReconstructionImageFilter.cxx
 * Sample DTI datasets may be obtained from 
 \begin verbatim
     ftp://public.kitware.com/pub/namic/DTI/Data/dwi.nhdr
     ftp://public.kitware.com/pub/namic/DTI/Data/dwi.img.gz ( gunzip this )
 \end verbatim
 *
 * \sa DiffusionTensor3D SymmetricSecondRankTensor
 * DiffusionTensor3DReconstructionLinearImageFilter
 * DiffusionTensor3DReconstructionNonlinearImageFilter
 * DiffusionTensor3DReconstructionWeightedImageFilter
 * DiffusionTensor3DReconstructionRicianImageFilter
 * \ingroup Multithreaded  TensorObjects
 */

template< class TGradientImagePixelType,
          class TTensorPrecision=double >
class ITK_EXPORT DiffusionTensor3DReconstructionImageFilterBase :
    public ImageToImageFilter< VectorImage< TGradientImagePixelType, 3 >, 
                               Image< DiffusionTensor3D< TTensorPrecision >, 3 > >
{
  
public:
  
  typedef DiffusionTensor3DReconstructionImageFilterBase Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef ImageToImageFilter< VectorImage< TGradientImagePixelType, 3 >, 
                              Image< DiffusionTensor3D< TTensorPrecision >, 3 > >
                                                  Superclass;
  
  /** Runtime information support. */
  itkTypeMacro(DiffusionTensor3DReconstructionImageFilterBase, 
                                                      ImageToImageFilter);
 
  typedef TGradientImagePixelType                  GradientPixelType;

  typedef typename Superclass::OutputImageType     OutputImageType;

  typedef typename OutputImageType::PixelType      TensorPixelType;

  typedef OutputImageType                          TensorImageType;

  typedef typename Superclass::OutputImageRegionType
                                                   OutputImageRegionType;

  /** Typedef defining one (of the many) gradient images.  */
  typedef Image< GradientPixelType, 3 >            ScalarImageType;

  /** An alternative typedef defining one (of the many) gradient images. 
   * It will be assumed that the vectorImage has the same dimension as the 
   * Gradient image and a vector length parameter of \c n (number of
   * gradient directions)*/
  typedef typename Superclass::InputImageType      GradientImagesType;

  /** Holds each magnetic field gradient used to acquire one DWImage */
  typedef vnl_vector_fixed< TTensorPrecision, 3 >  GradientDirectionType;

  /** Container to hold gradient directions of the 'n' DW measurements */
  typedef VectorContainer< unsigned int, 
          GradientDirectionType >                  GradientDirectionContainerType;
  

  /** Another set method to add a gradient directions and its corresponding
   * image. The image here is a VectorImage. The user is expected to pass the 
   * gradient directions in a container. The ith element of the container 
   * corresponds to the gradient direction of the ith component image the 
   * VectorImage.  For the baseline image, a vector of all zeros
   * should be set.*/
  virtual void SetGradientImage( GradientDirectionContainerType *, 
                                 const GradientImagesType *image);
  
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
  itkSetMacro( Threshold, GradientPixelType );
  itkGetMacro( Threshold, GradientPixelType );

  /** 
   * The BValue \f$ (s/mm^2) \f$ value used in normalizing the tensors to 
   * physically meaningful units.  See equation (24) of the first reference for
   * a description of how this is applied to the tensor estimation.
   * Equation (1) of the same reference describes the physical significance.
   */
  itkSetMacro( BValue, TTensorPrecision);
#ifdef GetBValue
#undef GetBValue
#endif
  itkGetConstReferenceMacro( BValue, TTensorPrecision);

  /** Set whether the baseline image should be estimated during the
  tensor estimation.  If set to true the baseline signal is available
  via the GetBaseline method.*/
  virtual void SetEstimateBaseline(bool eb);

  /** Returns the baseline signal if it was computed.  Otherwise
  returns NULL */
  virtual typename ScalarImageType::Pointer GetBaseline();
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(ReferenceEqualityComparableCheck,
    (Concept::EqualityComparable<GradientPixelType>));
  itkConceptMacro(TensorEqualityComparableCheck,
    (Concept::EqualityComparable<TensorPixelType>));
  itkConceptMacro(GradientConvertibleToDoubleCheck,
    (Concept::Convertible<GradientPixelType, TTensorPrecision>));
  itkConceptMacro(DoubleConvertibleToTensorCheck,
    (Concept::Convertible<TTensorPrecision, TensorPixelType>));
  itkConceptMacro(GradientReferenceAdditiveOperatorsCheck,
    (Concept::AdditiveOperators<GradientPixelType, GradientPixelType,
                                GradientPixelType>));
  itkConceptMacro(ReferenceOStreamWritableCheck,
    (Concept::OStreamWritable<GradientPixelType>));
  itkConceptMacro(TensorOStreamWritableCheck,
    (Concept::OStreamWritable<TensorPixelType>));
  /** End concept checking */
#endif

protected:
  DiffusionTensor3DReconstructionImageFilterBase();
  virtual ~DiffusionTensor3DReconstructionImageFilterBase() {};
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void ComputeTensorBasis();
  
  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData( const 
      OutputImageRegionType &outputRegionForThread, int);

  /** Derived classes should override this method to estimate the
  tensor from the diffusion weighted signal.  The gradient directions,
  b-matrix and inverse of b-matrix are available as class members.
  \param contains all the diffusion weighted signals casted into real
  type.
   \return vector which contain the tensor elements as the first 6
  elements in row-major order D_xx, D_xy, D_xz, D_yy, D_yz, D_zz.
  The seventh component is the estimated baseline b_0 signal.*/
  virtual vnl_vector< TTensorPrecision > 
    EstimateTensor(const vnl_vector<TTensorPrecision>& S) const = 0;
  

  /** Holds the tensor basis coefficients G_k */
  typedef vnl_matrix< TTensorPrecision >            TensorBasisMatrixType; 
  
  /** Matrix encoding of gradient directions and b-values.  This is
   * the matrix multiplied by the tensor and T2 signal which gives the
   * diffusion weighted signal. */
  TensorBasisMatrixType                             m_BMatrix;

  /** psuedo-inverse of m_BMatrix */
  TensorBasisMatrixType                             m_TensorBasis;
  
  /** container to hold gradient directions */
  typename GradientDirectionContainerType::Pointer  m_GradientDirectionContainer;

  /** Number of gradient measurements */
  unsigned int                                      m_NumberOfGradientDirections;

  /** Threshold on the reference image data */
  GradientPixelType                                 m_Threshold;

  /** LeBihan's b-value for normalizing tensors.  This is the b-value
   * for unit magnitude directions.  Multiple b-values are handled by
   * scaling of the gradient vectors. */
  TTensorPrecision                                  m_BValue;

private:
  /** Whether the baseline signal should be estimated and saved */
  bool                                              m_EstimateBaseline;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionImageFilterBase.txx"
#endif

#endif

