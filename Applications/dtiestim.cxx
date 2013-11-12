/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/03/03 15:15:31 $
  Version:   $Revision: 1.10 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// This program estimates a single diffusion tensor model at every
// voxel in an image.

// STL includes
#include <string>
#include <iostream>
#include <iomanip>

// ITK includes
// datastructures
#include <itkImage.h>
#include <itkVector.h>
#include <itkMetaDataObject.h>
#include <itkVersion.h>

// IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Filters
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkLogImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkExpImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkCastImageFilter.h>
#include <itkNthElementImageAdaptor.h>
#if ITK_VERSION_MAJOR < 4
#include <itkOtsuThresholdImageCalculator.h>
#else
#include <itkOtsuThresholdImageFilter.h>
#endif

#include "itkVectorMaskNegatedImageFilter.h"
#include "itkVectorMaskImageFilter.h"
#include "itkDiffusionTensor3DReconstructionNonlinearImageFilter.h"
#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.h"
#include "itkDiffusionTensor3DReconstructionRicianImageFilter.h"
#include "itkDiffusionTensor3DReconstructionLinearImageFilter.h"
#include "itkTensorRotateImageFilter.h"

#include <vnl/algo/vnl_svd.h>

#include "dtitypes.h"
#include <dtiestimCLP.h>

const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

enum EstimationType { LinearEstimate, NonlinearEstimate, WeightedEstimate, MaximumLikelihoodEstimate };

#if 0
void validate(boost::any& v,
              const std::vector<std::string>& values,
              EstimationType* target_type,
              int)
{
  using namespace boost::program_options;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if( s == "lls" || s == "linear" )
    {
    v = any(LinearEstimate);
    }
  else if( s == "nls" || s == "nonlinear" )
    {
    v = any(NonlinearEstimate);
    }
  else if( s == "wls" || s == "weighted" )
    {
    v = any(WeightedEstimate);
    }
  else if( s == "ml" )
    {
    v = any(MaximumLikelihoodEstimate);
    }
  else
    {
    throw validation_error("Estimation type invalid.  Only \"lls\", \"nls\", \"wls\", and \"ml\" allowed.");
    }
}

#endif

int main(int argc, char* argv[])
{
  PARSE_ARGS;
  // End option reading configuration

  // Display help if asked or program improperly called
  if( dwiImage == "" || tensorOutput == "" )
    {
    /*   if(help == true)
    {
      std::cout << "Version: $Date: 2009-03-03 15:15:31 $ $Revision: 1.10 $" << std::endl;
      std::cout << ITK_SOURCE_VERSION << std::endl;
      return EXIT_SUCCESS;
    }
    else
    {*/
    std::cerr << "DWI image and output tensor filename needs to be specified." << std::endl;
    return EXIT_FAILURE;
    /* }
    */
    }

  bool VERBOSE(verbose);
  if( stepSize < 0.0 )
    {
//   if(vm["method"].as<EstimationType>() == NonlinearEstimate ||
//      vm["method"].as<EstimationType>() == MaximumLikelihoodEstimate)
    if( method == "nls" || method == "ml" )
      {
      std::cerr << "Step size not set for optimization method" << std::endl;
      return EXIT_FAILURE;
      }
    }
  if( sigma < 0.0 )
    {
    //    if(vm["method"].as<EstimationType>() == MaximumLikelihoodEstimate)
    if( method == "ml" )
      {
      std::cerr << "Noise level not set for optimization method" << std::endl;
      return EXIT_FAILURE;
      }
    }
  // Read diffusion weighted MR

  typedef itk::ImageFileReader<VectorImageType> FileReaderType;
  FileReaderType::Pointer dwireader = FileReaderType::New();
  dwireader->SetFileName(dwiImage.c_str() );

  try
    {
    if( VERBOSE )
      {
      std::cout << "Reading Data" << std::endl;
      }
    dwireader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  VectorImageType::Pointer dwi = dwireader->GetOutput();

  // Read dwi meta-data

  // read into b0 the DWMRI_b-value that is the b-value of
  // the experiment
  double b0 = 0;
  bool   readbvalue = false;

  // read into gradientContainer the gradients
  typedef itk::DiffusionTensor3DReconstructionLinearImageFilter<ScalarPixelType, RealType>
    DiffusionEstimationFilterType;
  typedef itk::DiffusionTensor3DReconstructionNonlinearImageFilter<ScalarPixelType, RealType>
    NLDiffusionEstimationFilterType;
  typedef itk::DiffusionTensor3DReconstructionRicianImageFilter<ScalarPixelType, ScalarPixelType, RealType>
    MLDiffusionEstimationFilterType;
  typedef itk::DiffusionTensor3DReconstructionWeightedImageFilter<ScalarPixelType, RealType>
    WLDiffusionEstimationFilterType;

  DiffusionEstimationFilterType::GradientDirectionContainerType::Pointer gradientContainer =
    DiffusionEstimationFilterType::GradientDirectionContainerType::New();

  typedef DiffusionEstimationFilterType::GradientDirectionType GradientType;

  itk::MetaDataDictionary & dict = dwi->GetMetaDataDictionary();

  vnl_matrix<double> transform(3, 3);
  transform.set_identity();

  // Apply measurement frame if it exists
  if( dict.HasKey(NRRD_MEASUREMENT_KEY) )
    {
    if( VERBOSE )
      {
      std::cout << "Reorienting gradient directions to image coordinate frame" << std::endl;
      }

    // measurement frame
    vnl_matrix<double> mf(3, 3);
    // imaging frame
    vnl_matrix<double>                imgf(3, 3);
    std::vector<std::vector<double> > nrrdmf;
    itk::ExposeMetaData<std::vector<std::vector<double> > >(dict, NRRD_MEASUREMENT_KEY, nrrdmf);

    imgf = dwi->GetDirection().GetVnlMatrix();
    if( VERBOSE )
      {
      std::cout << "Image frame: " << std::endl;
      std::cout << imgf << std::endl;
      }
    for( unsigned int i = 0; i < 3; ++i )
      {
      for( unsigned int j = 0; j < 3; ++j )
        {
        mf(i, j) = nrrdmf[j][i];

        nrrdmf[j][i] = imgf(i, j);
        }
      }

    if( VERBOSE )
      {
      std::cout << "Meausurement frame: " << std::endl;
      std::cout << mf << std::endl;
      }

    itk::EncapsulateMetaData<std::vector<std::vector<double> > >(dict, NRRD_MEASUREMENT_KEY, nrrdmf);
    // prevent slicer error

    transform = vnl_svd<double>(imgf).inverse() * mf;

    if( VERBOSE )
      {
      std::cout << "Transform: " << std::endl;
      std::cout << transform << std::endl;
      }

    }

  if( dict.HasKey("modality") )
    {
    itk::EncapsulateMetaData<std::string>(dict, "modality", "DTMRI");
    }

  std::vector<std::string> keys = dict.GetKeys();
  for( std::vector<std::string>::const_iterator it = keys.begin();
       it != keys.end(); ++it )
    {
    if( it->find("DWMRI_b-value") != std::string::npos )
      {
      std::string t;
      itk::ExposeMetaData<std::string>(dict, *it, t);
      readbvalue = true;
      b0 = atof(t.c_str() );
      }
    else if( it->find("DWMRI_gradient") != std::string::npos )
      {
      std::string value;

      itk::ExposeMetaData<std::string>(dict, *it, value);
      std::istringstream iss(value);
      GradientType       g;
      iss >> g[0] >> g[1] >> g[2];

      g = transform * g;

      unsigned int ind;
      std::string  temp = it->substr(it->find_last_of('_') + 1);
      ind = atoi(temp.c_str() );

      gradientContainer->InsertElement(ind, g);
      }
    }
  for( std::vector<std::string>::const_iterator it = keys.begin();
       it != keys.end(); ++it )
    {
    if( it->find("DWMRI_NEX") != std::string::npos )
      {
      std::string numrepstr;

      itk::ExposeMetaData<std::string>(dict, *it, numrepstr);
      unsigned int numreps = atoi(numrepstr.c_str() );

      std::string  indtorepstr = it->substr(it->find_last_of('_') + 1);
      unsigned int indtorep =  atoi(indtorepstr.c_str() );

      GradientType g = gradientContainer->GetElement(indtorep);
      for( unsigned int i = indtorep + 1; i < indtorep + numreps; i++ )
        {
        gradientContainer->InsertElement(i, g);
        }
      }

    }

  if( VERBOSE )
    {
    std::cout << "NGrads: " << gradientContainer->Size() << std::endl;
    for( unsigned int i = 0; i < gradientContainer->Size(); ++i )
      {
      std::cout << gradientContainer->GetElement(i) << std::endl;
      }
    }

  if( !readbvalue )
    {
    std::cerr << "BValue not specified in header file" << std::endl;
    return EXIT_FAILURE;
    }

  if( VERBOSE )
    {
    std::cout << "BValue: " << b0 << std::endl;
    }

  if( VERBOSE )
    {
    std::cout << "Effective b-value per-direction" << std::endl;
    for( unsigned int i = 0; i < gradientContainer->Size(); ++i )
      {
      const double gsqnorm = gradientContainer->GetElement(i).squared_magnitude();
      std::cout << "G[" << std::setw(4) << std::setfill('0') << i << "]: " << b0 * gsqnorm << std::endl;
      }

    }

  // Read brain mask if it is specified.
  if( brainMask != "" )
    {
    typedef itk::ImageFileReader<LabelImageType> MaskFileReaderType;
    MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();
    maskreader->SetFileName(brainMask.c_str() );

    try
      {
      maskreader->Update();

      if( VERBOSE )
        {
        std::cout << "Masking Data" << std::endl;
        }

      // Check for same pixelsize & origin, and correct if necessary.
      LabelImageType::RegionType maskRegion = (maskreader->GetOutput())->GetLargestPossibleRegion();
      VectorImageType::RegionType dwiRegion = (dwireader->GetOutput())->GetLargestPossibleRegion();
      LabelImageType::SpacingType maskSpacing = (maskreader->GetOutput())->GetSpacing();
      VectorImageType::SpacingType dwiSpacing = (dwireader->GetOutput())->GetSpacing();
      LabelImageType::PointType maskOrigin = (maskreader->GetOutput())->GetOrigin();
      VectorImageType::PointType dwiOrigin = (dwireader->GetOutput())->GetOrigin();

      if ( ( (maskRegion.GetSize(0) == dwiRegion.GetSize(0)) && (maskRegion.GetSize(1) == dwiRegion.GetSize(1))
	     && (maskRegion.GetSize(2) == dwiRegion.GetSize(2)))) {
	if ( (maskSpacing[0] != dwiSpacing[0]) || (maskSpacing[1] != dwiSpacing[1]) || (maskSpacing[2] != dwiSpacing[2])) {
	  // Same size, but different pixeldim
	  std::cout << "warning pixelsizes are off between mask and dwi, ignoring and masking nevertheless" << std::endl;
	  (maskreader->GetOutput())->SetSpacing(dwiSpacing);
	}
	if ( (maskOrigin[0] != dwiOrigin[0]) || (maskOrigin[1] != dwiOrigin[1]) || (maskOrigin[2] != dwiOrigin[2])) {
	  // Same size, but different pixeldim
	  std::cout << "warning origin locations are off between mask and dwi, ignoring and masking nevertheless" << std::endl;
	  (maskreader->GetOutput())->SetOrigin(dwiOrigin);
	}
      }

      typedef itk::VectorMaskImageFilter<VectorImageType, LabelImageType, VectorImageType> MaskFilterType;
      MaskFilterType::Pointer mask = MaskFilterType::New();
//      mask->ReleaseDataFlagOn();
      mask->SetInput1(dwireader->GetOutput() );
      mask->SetInput2(maskreader->GetOutput() );
      mask->SetCoordinateTolerance( 0.001 ) ;
      mask->SetDirectionTolerance( 0.001 ) ;
      mask->Update();

      dwi = mask->GetOutput();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }

    }

  // Read negative mask
  if( badRegionMask != "" )
    {
    typedef itk::ImageFileReader<LabelImageType> MaskFileReaderType;
    MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();
    maskreader->SetFileName(badRegionMask);

    //  Go ahead and read data so we can use adaptors as necessary
    try
      {
      if( VERBOSE )
        {
        std::cout << "Masking Bad Regions" << std::endl;
        }

      maskreader->Update();

      typedef itk::VectorMaskNegatedImageFilter<VectorImageType, LabelImageType, VectorImageType> MaskFilterType;
      MaskFilterType::Pointer mask = MaskFilterType::New();
//      mask->ReleaseDataFlagOn();
      mask->SetInput1(dwi);
      mask->SetInput2(maskreader->GetOutput() );
      mask->Update();

      dwi = mask->GetOutput();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Compute average B0 image
  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, RealImageType> VectorSelectionFilterType;
  VectorSelectionFilterType::Pointer biextract = VectorSelectionFilterType::New();
  biextract->SetInput(dwi);
  int numberB0Directions = 0;

  RealImageType::Pointer B0Image;
  for( unsigned int directionIndex = 0; directionIndex < gradientContainer->Size(); directionIndex++ )
    {
    // get information whether image is b0 or not
    GradientType g = gradientContainer->GetElement(directionIndex);
    if( g[0] == 0 && g[1] == 0 && g[2] == 0 )
      {
      // image is not b0 image
      biextract->SetIndex(directionIndex);

      if( numberB0Directions == 0 )
        {
        B0Image = biextract->GetOutput();
        }
      else
        {
        // add it to the current B0 image
        try
          {
          typedef itk::AddImageFilter<RealImageType> AddImageFilterType;
          AddImageFilterType::Pointer addfilter = AddImageFilterType::New();
          addfilter->SetInput1(biextract->GetOutput() );
          addfilter->SetInput2(B0Image);
          addfilter->Update();
          B0Image = addfilter->GetOutput();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Error in addition for B0 computation" << std::endl;
          std::cerr << e << std::endl;
          }

        }

      numberB0Directions++;
      }
    }

  if( numberB0Directions == 0 )
    {
    if( VERBOSE )
      {
      std::cout << "No B0 image, setting first gradient image as B0 image (rather random behavior though)" << std::endl;
      }
    biextract->SetIndex(0);
    B0Image = biextract->GetOutput();
    }
  else
    {
    typedef itk::ImageRegionIterator<RealImageType> RealIterator;
    RealIterator iterImage(B0Image, B0Image->GetBufferedRegion() );
    while( !iterImage.IsAtEnd() )
      {
      iterImage.Set(iterImage.Get() / numberB0Directions);
      ++iterImage;
      }
    }

  if( VERBOSE )
    {
    std::cout << "Number of  B0 images : " << numberB0Directions << std::endl;
    }

  // Output B0 image if requested
  if( B0 != "" )
    {
    try
      {
      if( VERBOSE )
        {
        std::cout << "Writing B0" << std::endl;
        }
      typedef itk::ImageFileWriter<RealImageType> RealImageFileWriterType;
      RealImageFileWriterType::Pointer realwriter = RealImageFileWriterType::New();
      realwriter->SetInput(B0Image);
      realwriter->SetUseCompression(true);
      realwriter->SetFileName(B0.c_str() );
      realwriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Could not write B0 file" << std::endl;
      std::cerr << e << std::endl;
      }
    }

  // If we didnt specify a threshold compute it as the ostu threshold
  // of the baseline image

  ScalarPixelType _threshold;
  // the ScalarPixelType is unsigned short. ModuleDescription only has
  // 'int' as a possible type, which should be wider than unsigned short, and
  //  contains its range, which is 0..ScalarPixelType.max()
  //  So it's 'safe' to initialize threshold to -1 and use that as a sentinel
  //  value indicating that no threshold was specified on the command line.
  if( threshold >= 0 )
    {
    _threshold = static_cast<ScalarPixelType>(threshold);
    }
  else
    {
#if ITK_VERSION_MAJOR < 4
    typedef itk::OtsuThresholdImageCalculator<RealImageType>
      OtsuThresholdCalculatorType;

    OtsuThresholdCalculatorType::Pointer otsucalculator = OtsuThresholdCalculatorType::New();
    otsucalculator->SetImage(B0Image);
    otsucalculator->Compute();
#else
    typedef itk::OtsuThresholdImageFilter<RealImageType,RealImageType>
      OtsuThresholdImageFilterType;

    OtsuThresholdImageFilterType::Pointer otsucalculator = OtsuThresholdImageFilterType::New();
    otsucalculator->SetInput(B0Image);
    otsucalculator->Update();
#endif
    _threshold = static_cast<ScalarPixelType>(.9 * otsucalculator->GetThreshold() );
    if( VERBOSE )
      {
      std::cout << "Otsu threshold: " << threshold << std::endl;
      }
    }

  // Output b0 threshold mask if requested
  // BUG in original -- looked for "threshold-mask" tag in command line, when
  // it was really named B0
  if( B0MaskOutput != "" )
    {
    // Will take last B0 image in sequence

    typedef itk::BinaryThresholdImageFilter<RealImageType, LabelImageType> ThresholdFilterType;
    ThresholdFilterType::Pointer thresholdfilter = ThresholdFilterType::New();
    thresholdfilter->SetInput(B0Image);
    thresholdfilter->SetLowerThreshold(_threshold);
    thresholdfilter->SetUpperThreshold(itk::NumericTraits<ScalarPixelType>::max() );
    thresholdfilter->Update();

    try
      {
      if( VERBOSE )
        {
        std::cout << "Writing mask B0" << std::endl;
        }
      typedef itk::ImageFileWriter<LabelImageType> MaskImageFileWriterType;
      MaskImageFileWriterType::Pointer maskwriter = MaskImageFileWriterType::New();
      maskwriter->SetInput(thresholdfilter->GetOutput() );
      maskwriter->SetFileName(B0MaskOutput.c_str() );
      maskwriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Could not write threshold mask file" << std::endl;
      std::cerr << e << std::endl;
      }

    }

  // Output idwi image if requested
  if( IDWI != "" )
    {
    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, RealImageType> VectorSelectionFilterType;
    VectorSelectionFilterType::Pointer _biextract = VectorSelectionFilterType::New();
    _biextract->SetInput(dwi);
    int numberNonB0Directions = 0;

    RealImageType::Pointer idwiImage;

    typedef itk::LogImageFilter<RealImageType, RealImageType> LogImageFilterType;
    for( unsigned int directionIndex = 0; directionIndex < gradientContainer->Size(); directionIndex++ )
      {
      // get information whether image is b0 or not
      GradientType g = gradientContainer->GetElement(directionIndex);
      if( g[0] != 0 || g[1] != 0 || g[2] != 0 )
        {
        // image is not b0 image
        _biextract->SetIndex(directionIndex);

        if( numberNonB0Directions == 0 )
          {
          // log of the first image and set it as current idwi image
          try
            {
            LogImageFilterType::Pointer logfilter = LogImageFilterType::New();
            logfilter->SetInput(_biextract->GetOutput() );
            logfilter->Update();
            idwiImage = logfilter->GetOutput();
            }
          catch( itk::ExceptionObject & e )
            {
            std::cerr << "Error in log computation" << std::endl;
            std::cerr << e << std::endl;
            }

          }
        else
          {
          // log of the image and add it to the current idwi image
          try
            {
            LogImageFilterType::Pointer logfilter = LogImageFilterType::New();
            logfilter->SetInput(_biextract->GetOutput() );
            logfilter->Update();
            typedef itk::AddImageFilter<RealImageType> AddImageFilterType;
            AddImageFilterType::Pointer addfilter = AddImageFilterType::New();
            addfilter->SetInput1(logfilter->GetOutput() );
            addfilter->SetInput2(idwiImage);
            addfilter->Update();
            idwiImage = addfilter->GetOutput();
            }
          catch( itk::ExceptionObject & e )
            {
            std::cerr << "Error in log computation" << std::endl;
            std::cerr << e << std::endl;
            }

          }

        numberNonB0Directions++;
        }
      }

    // idwiImage contains the sum of all log transformed directional images
    // need to divide by numberNonB0Directions and compute exponential image

    if( VERBOSE )
      {
      std::cout << "Number of non B0 images : " << numberNonB0Directions << std::endl;
      }

    typedef itk::ImageRegionIterator<RealImageType> RealIterator;
    RealIterator iterImage(idwiImage, idwiImage->GetBufferedRegion() );
    while( !iterImage.IsAtEnd() )
      {
      iterImage.Set(iterImage.Get() / numberNonB0Directions);
      ++iterImage;
      }

    typedef itk::ExpImageFilter<RealImageType, RealImageType> ExpFilterType;
    ExpFilterType::Pointer expfilter = ExpFilterType::New();
    expfilter->SetInput(idwiImage);

    try
      {
      typedef itk::ImageFileWriter<RealImageType> RealImageFileWriterType;
      RealImageFileWriterType::Pointer realwriter = RealImageFileWriterType::New();
      realwriter->SetInput(expfilter->GetOutput() );
      realwriter->SetUseCompression(true);
      realwriter->SetFileName(IDWI.c_str() );
      realwriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Could not write idwi file" << std::endl;
      std::cerr << e << std::endl;
      }

    }

  // Estimate tensors
  typedef itk::ImageToImageFilter<VectorImageType, TensorImageType> DiffusionEstimationBaseType;
  TensorImageType::Pointer tensors;

  if( VERBOSE )
    {
    std::cout << "Estimation method: " << method << std::endl;
    }

  //  if(vm["method"].as<EstimationType>() == LinearEstimate)
  if( method == "lls" )
    {
       DiffusionEstimationFilterType::Pointer llsestimator = DiffusionEstimationFilterType::New();
       llsestimator->ReleaseDataFlagOn();
       llsestimator->SetGradientImage(gradientContainer, dwi);
       llsestimator->SetBValue(b0);
       llsestimator->SetThreshold(_threshold);
       llsestimator->Update();
       tensors = llsestimator->GetOutput();
    }
  //  else if(vm["method"].as<EstimationType>() == NonlinearEstimate)
  else if( method == "nls" )
    {
    NLDiffusionEstimationFilterType::Pointer estimator = NLDiffusionEstimationFilterType::New();
    estimator->ReleaseDataFlagOn();

    estimator->SetGradientImage(gradientContainer, dwi);
    estimator->SetBValue(b0);
    estimator->SetThreshold(_threshold);
    estimator->SetStep(stepSize);
    estimator->SetNumberOfThreads(1);
    estimator->Update();
    tensors = estimator->GetOutput();
    }
  // else if(vm["method"].as<EstimationType>() == WeightedEstimate)
  else if( method == "wls" )
    {
    WLDiffusionEstimationFilterType::Pointer estimator = WLDiffusionEstimationFilterType::New();
    estimator->ReleaseDataFlagOn();

    if( VERBOSE )
      {
      std::cout << "Weighting steps: " << weightIterations << std::endl;
      }

    estimator->SetGradientImage(gradientContainer, dwi);
    estimator->SetBValue(b0);
    estimator->SetThreshold(_threshold);
    estimator->SetNumberOfIterations(weightIterations);
    estimator->Update();
    tensors = estimator->GetOutput();
    }
  //  else if(vm["method"].as<EstimationType>() == MaximumLikelihoodEstimate)
  else if( method == "ml" )
    {
    WLDiffusionEstimationFilterType::Pointer estimatorInit = WLDiffusionEstimationFilterType::New();
    estimatorInit->ReleaseDataFlagOn();

    estimatorInit->SetGradientImage(gradientContainer, dwi);
    estimatorInit->SetBValue(b0);
    estimatorInit->SetThreshold(_threshold);
    estimatorInit->SetNumberOfIterations(weightIterations);
    estimatorInit->Update();
    TensorImageType::Pointer inittensors = estimatorInit->GetOutput();

    MLDiffusionEstimationFilterType::Pointer estimator = MLDiffusionEstimationFilterType::New();
    estimator->ReleaseDataFlagOn();

    if( VERBOSE )
      {
      std::cout << "Weighting steps: " << weightIterations << std::endl;
      }

    estimator->SetGradientImage(gradientContainer, dwi);
    estimator->SetBValue(b0);
    estimator->SetThreshold(_threshold);
    estimator->SetInitialTensor(inittensors);
    estimator->SetStep(stepSize);
    estimator->SetNumberOfThreads(1);
    std::cout << "Start sigma: " << sigma << std::endl;
    estimator->SetSigma(sigma);
    estimator->Update();
    tensors = estimator->GetOutput();
    }
  else
    {
    std::cerr << "Invalid estimation method"  << std::endl;
    return EXIT_FAILURE;
    }

  // wp = D*x
  // wv = M*x
  // wv' = D'M*x

  // Write tensor file if requested
  try
    {
    if( !doubleDTI )
      {
      typedef itk::DiffusionTensor3D<float> TensorFloatPixelType ;
      typedef itk::Image<TensorFloatPixelType, DIM> TensorFloatImageType ;
      typedef itk::CastImageFilter< TensorImageType, TensorFloatImageType > CastDTIFilterType ;
      CastDTIFilterType::Pointer castFilter = CastDTIFilterType::New() ;
      castFilter->SetInput( tensors ) ;
      typedef itk::ImageFileWriter<TensorFloatImageType> TensorFileWriterType ;
      TensorFileWriterType::Pointer tensorWriter = TensorFileWriterType::New() ;
      tensorWriter->SetFileName(tensorOutput.c_str()) ;
      tensors->SetMetaDataDictionary(dict) ;
      tensorWriter->SetInput(castFilter->GetOutput()) ;
      tensorWriter->SetUseCompression(true) ;
      tensorWriter->Update() ;
      }
    else
      {
      typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType ;
      TensorFileWriterType::Pointer tensorWriter = TensorFileWriterType::New() ;
      tensorWriter->SetFileName(tensorOutput.c_str()) ;
      tensors->SetMetaDataDictionary(dict) ;
      tensorWriter->SetInput(tensors) ;
      tensorWriter->SetUseCompression(true) ;
      tensorWriter->Update() ;
      }
    }
  catch( itk::ExceptionObject e )
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
