/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: TensorEstimateTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007-09-04 20:12:29 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkNewDiffusionTensor3DReconstructionImageFilter.h"
#include "itkDiffusionTensor3DReconstructionNonlinearImageFilter.h"
#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>

int TensorEstimateTest(int, char*[])
{
  typedef short int          DWIComponentPixelType;
  typedef double             TensorPrecisionType;

  typedef itk::NewDiffusionTensor3DReconstructionImageFilter< 
      DWIComponentPixelType, DWIComponentPixelType, TensorPrecisionType > 
        TensorReconstructionImageFilterType;
  typedef itk::DiffusionTensor3DReconstructionNonlinearImageFilter< 
      DWIComponentPixelType, DWIComponentPixelType, TensorPrecisionType > 
        NonlinearTensorReconstructionImageFilterType;
  typedef itk::DiffusionTensor3DReconstructionWeightedImageFilter< 
      DWIComponentPixelType, DWIComponentPixelType, TensorPrecisionType > 
        WeightedTensorReconstructionImageFilterType;


  typedef TensorReconstructionImageFilterType::GradientImagesType DWIImageType;

  // LLS
  TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = 
    TensorReconstructionImageFilterType::New();
  // NLS
  NonlinearTensorReconstructionImageFilterType::Pointer nlstensorReconstructionFilter = 
    NonlinearTensorReconstructionImageFilterType::New();
  // WLS
  WeightedTensorReconstructionImageFilterType::Pointer wlstensorReconstructionFilter = 
    WeightedTensorReconstructionImageFilterType::New();
  
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

  typedef itk::VectorContainer< unsigned int, vnl_vector_fixed<double, 3> > GradientContainer;
  GradientContainer::Pointer gradientDirections = GradientContainer::New();
  std::cout << "Gradients: " << std::endl;
  for(unsigned int i = 0; i < numberOfGradientImages ; ++i)
    {
    vnl_vector_fixed<double, 3> grad(PRIVATEGradientDirections[i]);
    gradientDirections->CastToSTLContainer().push_back(grad);
    std::cout << grad << std::endl;
    }
  
  // Create gradient images
  //
  typedef DWIImageType::Pointer GradientImagePointer;
  typedef DWIImageType::RegionType  GradientRegionType;
  typedef GradientRegionType::IndexType  GradientIndexType;
  typedef GradientRegionType::SizeType   GradientSizeType;
  typedef DWIImageType::PixelType   VectorDWIPixelType;
  
  DWIImageType::Pointer dwimage = DWIImageType::New();
  GradientSizeType  sizeGradientImage  = {{ 4, 4, 4 }};
  GradientIndexType indexGradientImage = {{ 0, 0, 0 }};
  GradientRegionType     regionGradientImage;
  regionGradientImage.SetSize(  sizeGradientImage );
  regionGradientImage.SetIndex( indexGradientImage);
  dwimage->SetVectorLength(numberOfGradientImages);
  dwimage->SetRegions( regionGradientImage );
  dwimage->Allocate();
  
  const double PRIVATEtruetensor[6] = {3.0e-4, 0.0, 0.0, 2.0e-4, 0.0, 1.0e-4};
  itk::DiffusionTensor3D<double> truetensor(PRIVATEtruetensor);

  itk::ImageRegionIteratorWithIndex< DWIImageType > git( 
    dwimage, regionGradientImage );
  git.GoToBegin();
  while( !git.IsAtEnd() )
    {
    VectorDWIPixelType dwiVector = git.Get();
    dwiVector[0] = 100; // Reference intensity of 100
    for(unsigned int i = 1; i < numberOfGradientImages; ++i)
      {
      std::cout << "ADC[" << i << "]: " << truetensor.GetApparentDiffusionCoefficient(gradientDirections->ElementAt(i)) << std::endl;
      std::cout << -bvalue*truetensor.GetApparentDiffusionCoefficient(gradientDirections->ElementAt(i)) << std::endl;
      std::cout << dwiVector[0]*exp(-bvalue*truetensor.GetApparentDiffusionCoefficient(gradientDirections->ElementAt(i))) << std::endl;
      dwiVector[i] = static_cast<short int>(dwiVector[0]*exp(-bvalue*truetensor.GetApparentDiffusionCoefficient(gradientDirections->ElementAt(i))));
      }
    git.Set( dwiVector );
    ++git;
    }
    
  TensorReconstructionImageFilterType::GradientDirectionType gradientDirection;
  tensorReconstructionFilter->SetGradientImage( gradientDirections, dwimage );   
  nlstensorReconstructionFilter->SetGradientImage( gradientDirections, dwimage );
  wlstensorReconstructionFilter->SetGradientImage( gradientDirections, dwimage );

  tensorReconstructionFilter->SetBValue(bvalue);
  tensorReconstructionFilter->Update();

  nlstensorReconstructionFilter->SetInitialTensor(tensorReconstructionFilter->GetOutput());
  nlstensorReconstructionFilter->SetBValue(bvalue);
  nlstensorReconstructionFilter->Update();

  wlstensorReconstructionFilter->SetInitialTensor(tensorReconstructionFilter->GetOutput());
  wlstensorReconstructionFilter->SetBValue(bvalue);
  wlstensorReconstructionFilter->Update();

  typedef TensorReconstructionImageFilterType::TensorImageType TensorImageType;
  TensorImageType::Pointer tensorImage = tensorReconstructionFilter->GetOutput();
  TensorImageType::Pointer nlstensorImage = nlstensorReconstructionFilter->GetOutput();
  TensorImageType::Pointer wlstensorImage = wlstensorReconstructionFilter->GetOutput();
  typedef TensorImageType::IndexType TensorImageIndexType;
  
  TensorImageIndexType tensorImageIndex    = {{3,3,3}};

  std::cout << std::endl << "Pixels at index: " << tensorImageIndex << std::endl;
  
  bool llspassed = true;
  bool nlspassed = true;
  bool wlspassed = true;

  std::cout << "Diffusion weigthed vector" << std::endl;
  std::cout << dwimage->GetPixel(tensorImageIndex) << std::endl;

  itk::DiffusionTensor3D<double> llsresult(tensorImage->GetPixel(tensorImageIndex));
  itk::DiffusionTensor3D<double> nlsresult(nlstensorImage->GetPixel(tensorImageIndex));
  itk::DiffusionTensor3D<double> wlsresult(wlstensorImage->GetPixel(tensorImageIndex));

  double precision = 0.0001;
  for(unsigned int i = 0; i < 6; ++i)
    {
    llspassed = (vnl_math_abs(llsresult[i] - truetensor[i]) < precision);
    nlspassed = (vnl_math_abs(nlsresult[i] - truetensor[i]) < precision);
    wlspassed = (vnl_math_abs(wlsresult[i] - truetensor[i]) < precision);
    }

  bool failed = false;
  std::cout << std::endl << "LLS Reconstructed tensor : " << std::endl; 
  std::cout << llsresult << std::endl;
  if( !llspassed ) 
    {
    std::cout << "[FAILED]" << std::endl;
    
    std::cout << "Expected tensor : " << std::endl;
    std::cout << truetensor << std::endl;
    failed = true;
    }
  else
    {
    std::cout << "[PASSED]" << std::endl;
    }

  std::cout << std::endl << "NLS Reconstructed tensor : " << std::endl; 
  std::cout << nlsresult << std::endl;
  if( !nlspassed ) 
    {
    std::cout << "[FAILED]" << std::endl;
    
    std::cout << "Expected tensor : " << std::endl;
    std::cout << truetensor << std::endl;
    failed = true;
    }
  else
    {
    std::cout << "[PASSED]" << std::endl;
    }


  std::cout << std::endl << "WLS Reconstructed tensor : " << std::endl; 
  std::cout << wlsresult << std::endl;
  if( !wlspassed ) 
    {
    std::cout << "[FAILED]" << std::endl;
    
    std::cout << "Expected tensor : " << std::endl;
    std::cout << truetensor << std::endl;
    failed = true;
    }
  else
    {
    std::cout << "[PASSED]" << std::endl;
    }
      
  if(failed)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}

