/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: TensorEstimateTest.cxx,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkDiffusionTensor3DReconstructionLinearImageFilter.h"
#include "itkDiffusionTensor3DReconstructionNonlinearImageFilter.h"
#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.h"
#include "itkTensorApparentDiffusionCoefficientImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>

namespace{
  typedef short unsigned int  DWIComponentPixelType;
  typedef double              TensorPrecisionType;
  typedef itk::DiffusionTensor3D<TensorPrecisionType> TensorPixelType;
  typedef itk::Image<TensorPixelType, 3> TensorImageType;
  typedef itk::Image<DWIComponentPixelType, 3> ScalarImageType;

  const unsigned int numberOfGradientImages = 22;
  // Assign gradient directions
  //
  const double bvalue = 1000;
  const double PRIVATEGradientDirections[numberOfGradientImages][3] = 
  {
    {0, 0, 0},
    {-0.99823, -0.036499, 0.047033},
    {0.037648, 0.85056, 0.52452},
    {-0.85783, 0.49589, 0.13495},
    {0.25729, -0.025629, 0.96599},
    {-0.59516, -0.43791, 0.67381},
    {-0.54816, 0.49097, 0.67711},
    {-0.82631, 0.04477, 0.56144},
    {0.48971, 0.84011, 0.23324},
    {-0.82302, -0.53713, 0.18476},
    {0.77956, -0.037959, 0.62517},
    {-0.45054, 0.84917, 0.27554},
    {-0.3338, -0.024468, 0.94233},
    {-0.044461, 0.45427, 0.88975},
    {0.25947, -0.87224, 0.41457},
    {-0.33418, -0.85152, 0.40402},
    {-0.063334, -0.54422, 0.83655},
    {0.47387, -0.48088, 0.7377},
    {0.84672, 0.43078, 0.31224},
    {0.48855, 0.47013, 0.73505},
    {0.78015, -0.55755, 0.28372},
    {0.99647, -0.046849, 0.0696}
  };
};

template< class TTensorEstimateType >
class TensorEstimateTester
{
public:
  void initialize( typename TTensorEstimateType::Pointer testim)
  {
    typedef itk::VectorContainer< unsigned int, vnl_vector_fixed<double, 3> > GradientContainer;
    typename GradientContainer::Pointer gradientDirections = GradientContainer::New();
    std::cout << "Gradients: " << std::endl;
    for(unsigned int i = 0; i < numberOfGradientImages ; ++i)
    {
      vnl_vector_fixed<double, 3> grad(PRIVATEGradientDirections[i]);
      gradientDirections->CastToSTLContainer().push_back(grad);
      std::cout << grad << std::endl;
    }
    
    // Create gradient images
    //
    typedef typename TTensorEstimateType::GradientImagesType DWIImageType;
    typedef typename DWIImageType::Pointer GradientImagePointer;
    typedef typename DWIImageType::RegionType  GradientRegionType;
    typedef typename GradientRegionType::IndexType  GradientIndexType;
    typedef typename GradientRegionType::SizeType   GradientSizeType;
    typedef typename DWIImageType::PixelType   VectorDWIPixelType;
    
    typename DWIImageType::Pointer dwimage = DWIImageType::New();
    GradientSizeType  sizeGradientImage  = {{ 4, 4, 4 }};
    GradientIndexType indexGradientImage = {{ 0, 0, 0 }};
    GradientRegionType     regionGradientImage;
    regionGradientImage.SetSize(  sizeGradientImage );
    regionGradientImage.SetIndex( indexGradientImage);
    dwimage->SetVectorLength(numberOfGradientImages);
    dwimage->SetRegions( regionGradientImage );
    dwimage->Allocate();
  
    const double PRIVATEtruetensor[6] = {3.0e-4, 0.0, 0.0, 2.0e-4, 0.0, 1.0e-4};
    m_TrueTensor = TensorPixelType(PRIVATEtruetensor);
    
    itk::ImageRegionIteratorWithIndex< DWIImageType > git( 
      dwimage, regionGradientImage );
    git.GoToBegin();
    
    typedef itk::Functor::TensorApparentDiffusionCoefficient<itk::DiffusionTensor3D<double>, vnl_vector_fixed<double, 3>, double> ADC;
    while( !git.IsAtEnd() )
    {
      VectorDWIPixelType dwiVector = git.Get();
      m_B0 = dwiVector[0] = 1000; // Reference intensity of 100
      for(unsigned int i = 1; i < numberOfGradientImages; ++i)
      {
        dwiVector[i] = static_cast<short int>(dwiVector[0]*exp(-bvalue*ADC()(m_TrueTensor, gradientDirections->ElementAt(i))));
      }
      git.Set( dwiVector );
      ++git;
    }

    std::cout << "DWI signal: " << std::endl;
    typename DWIImageType::IndexType ind = {{3,3,3}};
    std::cout << dwimage->GetPixel(ind) << std::endl;
    
    testim->SetGradientImage( gradientDirections, dwimage );
    testim->SetBValue(bvalue);
    testim->Update();
  }

  bool checkB0Estimate( ScalarImageType::Pointer b0)
  {
    typedef TensorImageType::IndexType TensorImageIndexType;
  
    TensorImageIndexType tensorImageIndex    = {{3,3,3}};
    DWIComponentPixelType result(b0->GetPixel(tensorImageIndex));

    std::cout << "Estimated B0: " << result << std::endl;
    if (result != m_B0)
    {
        std::cout << "[FAILED]" << std::endl;
        std::cout << "Expected b0 : " << m_B0 << std::endl;
        return false;
    }
    std::cout << "[PASSED]" << std::endl;
    return true;
  }

  bool checkTensorEstimate( TensorImageType::Pointer timg)
  {
    typedef TensorImageType::IndexType TensorImageIndexType;
  
    TensorImageIndexType tensorImageIndex    = {{3,3,3}};

    std::cout << timg->GetPixel(tensorImageIndex) << std::endl;
  
    itk::DiffusionTensor3D<double> result(timg->GetPixel(tensorImageIndex));
    const double precision = 0.0001;
    for(unsigned int i = 0; i < 6; ++i)
    {
      if(fabs(result[i] - m_TrueTensor[i]) > precision)
      {
        std::cout << "[FAILED]" << std::endl;
        std::cout << "Expected tensor : " << std::endl;
        std::cout << m_TrueTensor << std::endl;
        return false;
      }
    }
    std::cout << "[PASSED]" << std::endl;
    return true;
  }

private:
  itk::DiffusionTensor3D<double> m_TrueTensor;
  DWIComponentPixelType m_B0;
};

int NLSTensorEstimateTest(int, char*[])
{
  typedef itk::DiffusionTensor3DReconstructionNonlinearImageFilter< 
        DWIComponentPixelType, TensorPrecisionType > 
        TensorReconstructionImageFilterType;

  typedef TensorReconstructionImageFilterType TensorReconstructionImageFilterType;

  // NLS
  TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = 
    TensorReconstructionImageFilterType::New();

  TensorEstimateTester<TensorReconstructionImageFilterType> tensortester;
  std::cout << "NLS Estimate" << std::endl;
  tensortester.initialize(tensorReconstructionFilter);

  std::cout << "NLS Reconstructed tensor" << std::endl;
  bool passed = tensortester.checkTensorEstimate(tensorReconstructionFilter->GetOutput());
  if(passed)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

int WLSTensorEstimateTest(int, char*[])
{
  typedef itk::DiffusionTensor3DReconstructionWeightedImageFilter< 
      DWIComponentPixelType, TensorPrecisionType > 
        WeightedTensorReconstructionImageFilterType;

  typedef WeightedTensorReconstructionImageFilterType TensorReconstructionImageFilterType;

  // WLS
  WeightedTensorReconstructionImageFilterType::Pointer wlstensorReconstructionFilter = 
    WeightedTensorReconstructionImageFilterType::New();

  TensorEstimateTester<WeightedTensorReconstructionImageFilterType> tensortester;
  std::cout << "WLS Estimate" << std::endl;
  tensortester.initialize(wlstensorReconstructionFilter);

  std::cout << "WLS Reconstructed tensor" << std::endl;
  bool passed = tensortester.checkTensorEstimate(wlstensorReconstructionFilter->GetOutput());
  if(passed)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

int LLSTensorEstimateTest(int, char*[])
{
  typedef itk::DiffusionTensor3DReconstructionLinearImageFilter< 
      DWIComponentPixelType, TensorPrecisionType > 
        TensorReconstructionImageFilterType;

  typedef TensorReconstructionImageFilterType TensorReconstructionImageFilterType;

  // WLS
  TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = 
    TensorReconstructionImageFilterType::New();

  TensorEstimateTester<TensorReconstructionImageFilterType> tensortester;
  std::cout << "LLS Estimate" << std::endl;
  tensorReconstructionFilter->SetEstimateBaseline(true);
  tensortester.initialize(tensorReconstructionFilter);

  std::cout << "LLS Reconstructed tensor" << std::endl;
  bool passed = tensortester.checkTensorEstimate(tensorReconstructionFilter->GetOutput());
  passed == passed && tensortester.checkB0Estimate(tensorReconstructionFilter->GetBaseline());
  if(passed)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

