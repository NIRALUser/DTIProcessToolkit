// This program calculates

// STL includes
#include <string>
#include <iostream>

// boost includes
#include <boost/program_options/option.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/cmdline.hpp>

// ITK includes
// datastructures
#include <itkImage.h>
#include <itkVector.h>
#include <itkMetaDataObject.h>

// IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Filters
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

#include <itkNthElementImageAdaptor.h>
#include <itkOtsuThresholdImageCalculator.h>

#include "itkVectorMaskNegatedImageFilter.h"
#include "itkVectorMaskImageFilter.h"
#include "itkDiffusionTensor3DReconstructionNonlinearImageFilter.h"
#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.h"
#include "itkDiffusionTensor3DReconstructionRicianImageFilter.h"
#include "itkNewDiffusionTensor3DReconstructionImageFilter.h"
#include "itkTensorRotateImageFilter.h"

#include <vnl/algo/vnl_svd.h>

const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

namespace po = boost::program_options;

enum EstimationType {Linear, Nonlinear, Weighted, MaximumLikelihood};
void validate(boost::any& v,
              const std::vector<std::string>& values,
              EstimationType* target_type,
              int)
{
  using namespace po;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s == "lls" || s == "linear")
    {
    v = any(Linear);
    }
  else if (s == "nls" || s == "nonlinear")
    {
    v = any(Nonlinear);
    }
  else if (s == "wls" || s == "weighted")
    {
    v = any(Weighted);
    }
  else if (s == "ml")
    {
    v = any(MaximumLikelihood);
    }
    else
    {
    throw validation_error("Estimation type invalid.  Only \"lls\", \"nls\", \"wls\", and \"ml\" allowed.");
    }
}

int main(int argc, char* argv[])
{
  // Define necessary types for images
  typedef double RealType;
  const unsigned int DIM = 3;
  typedef unsigned short DWIPixelType;
  typedef itk::VectorImage<DWIPixelType, DIM> VectorImageType;
  typedef itk::Image<unsigned char, DIM> MaskImageType;
  typedef itk::Image<RealType, DIM> RealImageType;
  typedef itk::Image<DWIPixelType, DIM> IntImageType;
  typedef itk::Image<itk::DiffusionTensor3D<double>, DIM> TensorImageType;

//  unsigned int scale;

  // Read program options/configuration
  po::options_description config("Usage: dtiestim dwi-image tensor-output [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("brain-mask,M", po::value<std::string>(), "sets brain mask")
    ("bad-region-mask,B", po::value<std::string>(), "sets bad region mask")
    ("threshold-mask,T", po::value<std::string>(), "File to write mask estimated from b0 and threshold")
    
     //("double,d", "Writes tensors in double precision -- currently default")
    ("threshold,t", po::value<DWIPixelType>(),"Baseline threshold for estimation")
    ("method,m", po::value<EstimationType>()->default_value(Linear,"Linear Method"),"Estimation method")
    ("step,s", po::value<double>()->default_value(1.0e-8),"Gradient descent step size")
    ("sigma", po::value<double>(),"Sigma for Rician ML estimation")
    ("verbose,v", "Verbose output")
  ;
//    ("fa-scale,s", po::value<unsigned int>(&scale)->default_value(10000),"FA scale factor.  If set the FA value is scaled by this factor and written out in an integer image format.")

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("dwi-image", po::value<std::string>(), "DWI image volume")
    ("tensor-output", po::value<std::string>(), "Tensor output")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("dwi-image",1);
  p.add("tensor-output",1);

  po::variables_map vm;

  try
    {
    po::store(po::command_line_parser(argc, argv).
            options(all).positional(p).run(), vm);
    po::notify(vm);     
    } 
  catch (const po::error &e)
    {
    std::cout << config << std::endl;
    return EXIT_FAILURE;
    }

  // End option reading configuration

  // Display help if asked or program improperly called
  if(vm.count("help") || !vm.count("dwi-image") || !vm.count("tensor-output"))
    {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      {
      std::cerr << "DWI image and output tensor filename needs to be specified." << std::endl;
      return EXIT_FAILURE;
      }
    }

  bool VERBOSE = false;
  if(vm.count("verbose"))
    VERBOSE = true;

  double step = 1.0e-8, sigma = 0.0;
  try
    {
    step = vm["step"].as<double>();
    }
  catch( ... )
    {
    if(vm["method"].as<EstimationType>() == Nonlinear || vm["method"].as<EstimationType>() == MaximumLikelihood)
      {
      std::cerr << "Step size not set for optimization method" << std::endl;
      return EXIT_FAILURE;
      }    
    }
  try
    {
    sigma = vm["sigma"].as<double>();
    }
  catch( ... )
    {
    if(VERBOSE)
      std::cout << "Sigma not set" << std::endl;
    if(vm["method"].as<EstimationType>() == MaximumLikelihood)
      {
      std::cerr << "Noise level not set for optimization method" << std::endl;
      return EXIT_FAILURE;
      }    
    
    }

  // Read diffusion weighted MR

  typedef itk::ImageFileReader<VectorImageType> FileReaderType;
  FileReaderType::Pointer dwireader = FileReaderType::New();
  dwireader->SetFileName(vm["dwi-image"].as<std::string>().c_str());

  try
    {
    dwireader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e <<std::endl;
    return EXIT_FAILURE;
    }

  VectorImageType::Pointer dwi = dwireader->GetOutput();

  // Read dwi meta-data

  // read into b0 the DWMRI_b-value that is the b-value of
  // the experiment
  double b0 = 0;
  bool readb0 = false;

  // read into gradientContainer the gradients
  typedef itk::NewDiffusionTensor3DReconstructionImageFilter<DWIPixelType,DWIPixelType, RealType>
    DiffusionEstimationFilterType;
  typedef itk::DiffusionTensor3DReconstructionNonlinearImageFilter<DWIPixelType,DWIPixelType, RealType>
    NLDiffusionEstimationFilterType;
  typedef itk::DiffusionTensor3DReconstructionRicianImageFilter<DWIPixelType,DWIPixelType, RealType>
    MLDiffusionEstimationFilterType;
  typedef itk::DiffusionTensor3DReconstructionWeightedImageFilter<DWIPixelType,DWIPixelType, RealType>
    WLDiffusionEstimationFilterType;


  DiffusionEstimationFilterType::GradientDirectionContainerType::Pointer gradientContainer = 
    DiffusionEstimationFilterType::GradientDirectionContainerType::New();
  
  typedef DiffusionEstimationFilterType::GradientDirectionType GradientType;

  itk::MetaDataDictionary & dict = dwi->GetMetaDataDictionary();
  std::vector<std::string> keys = dict.GetKeys();
  for(std::vector<std::string>::const_iterator it = keys.begin();
      it != keys.end(); ++it)
    {
    std::string value;
    if( it->find("DWMRI_b-value") != std::string::npos)
      {
      std::string t;
      itk::ExposeMetaData<std::string>(dict, *it, t);
      readb0 = true;
      b0 = atof(t.c_str());
      }
    else if( it->find("DWMRI_gradient") != std::string::npos)
      {
      std::string value;

      itk::ExposeMetaData<std::string>(dict, *it, value);
      std::istringstream iss(value);
      GradientType g;
      iss >> g[0] >> g[1] >> g[2];

//      g = g/ sqrt(2.0);

      unsigned int ind;
      std::string temp = it->substr(it->find_last_of('_')+1);
      ind = atoi(temp.c_str());
      
      gradientContainer->InsertElement(ind,g);
      }
    else if( it->find("DWMRI_NEX") != std::string::npos)
      {
      std::string numrepstr;

      itk::ExposeMetaData<std::string>(dict, *it, numrepstr);
      unsigned int numreps = atoi(value.c_str());

      std::string indtorepstr = it->substr(it->find_last_of('_')+1);
      unsigned int indtorep =  atoi(indtorepstr.c_str());

      GradientType g = gradientContainer->GetElement(indtorep);

      for(unsigned int i = indtorep+1; i < indtorep+numreps-1; i++)
        gradientContainer->InsertElement(i,g);
      }

    }

  if(!readb0)
    {
    std::cerr << "BValue not specified in header file" << std::endl;
    return EXIT_FAILURE;
    }

  if(VERBOSE)
    std::cout << "BValue: " << b0 << std::endl;

 //  if(vm.count("brain-mask") && vm.count("bad-region-mask"))
//     {
//     std::cerr << "--brain-mask and --bad-region-mask are mutually exclusive" << std::endl;
//     return EXIT_FAILURE;
//     }

  // Read brain mask if it is specified.  
  if(vm.count("brain-mask"))
    {
    typedef itk::ImageFileReader<MaskImageType> MaskFileReaderType;
    MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();
    maskreader->SetFileName(vm["brain-mask"].as<std::string>().c_str());

    try
      {
      maskreader->Update();
      
      typedef itk::VectorMaskImageFilter<VectorImageType,MaskImageType,VectorImageType> MaskFilterType;
      MaskFilterType::Pointer mask = MaskFilterType::New();
      mask->SetInput1(dwireader->GetOutput());
      mask->SetInput2(maskreader->GetOutput());
      mask->Update();

      dwi = mask->GetOutput();
      }
    catch (itk::ExceptionObject & e)
      {
      std::cerr << e <<std::endl;
      return EXIT_FAILURE;
      }

    }
  
  // Read negative mask
  if(vm.count("bad-region-mask"))
    {
    typedef itk::ImageFileReader<MaskImageType> MaskFileReaderType;
    MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();
    maskreader->SetFileName(vm["bad-region-mask"].as<std::string>().c_str());
    
    //  Go ahead and read data so we can use adaptors as necessary
    try
      {
      maskreader->Update();
      
      typedef itk::VectorMaskNegatedImageFilter<VectorImageType,MaskImageType,VectorImageType> MaskFilterType;
      MaskFilterType::Pointer mask = MaskFilterType::New();
      mask->SetInput1(dwi);
      mask->SetInput2(maskreader->GetOutput());
      mask->Update();

      dwi = mask->GetOutput();
      }
    catch (itk::ExceptionObject & e)
      {
      std::cerr << e <<std::endl;
      return EXIT_FAILURE;
      }
    }

  // If we didnt specify a threshold compute it as the ostu threshold
  // of the baseline image
  

  DWIPixelType threshold;
  if(vm.count("threshold"))
    {
    threshold = vm["threshold"].as<DWIPixelType>();
    }
  else
    {
    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, IntImageType>
      BaselineExtractAdaptorType;

    BaselineExtractAdaptorType::Pointer baselineextract = BaselineExtractAdaptorType::New();
    baselineextract->SetInput(dwireader->GetOutput());
    baselineextract->SetIndex(0);
    baselineextract->Update();

    typedef itk::OtsuThresholdImageCalculator<IntImageType> 
      OtsuThresholdCalculatorType;

    OtsuThresholdCalculatorType::Pointer  otsucalculator = OtsuThresholdCalculatorType::New();
    otsucalculator->SetImage(baselineextract->GetOutput());
    otsucalculator->Compute();
    threshold = static_cast<DWIPixelType>(.9 * otsucalculator->GetThreshold());

    if(VERBOSE)
      std::cout << "Otsu threshold: " << threshold << std::endl;
    }

  // Output b0 threshold mask if requested
  if(vm.count("threshold-mask"))
    {
    // TODO: this will not work if the baselin image is not uniquely
    // the first image in the stack

    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType,IntImageType> VectorSelectionFilterType;
    VectorSelectionFilterType::Pointer b0extract = VectorSelectionFilterType::New();
    b0extract->SetInput(dwi);
    b0extract->SetIndex(0);

    typedef itk::BinaryThresholdImageFilter<IntImageType,MaskImageType> ThresholdFilterType;
    ThresholdFilterType::Pointer thresholdfilter = ThresholdFilterType::New();
    thresholdfilter->SetInput(b0extract->GetOutput());
    thresholdfilter->SetLowerThreshold(threshold);
    thresholdfilter->SetUpperThreshold(itk::NumericTraits<DWIPixelType>::max());
    thresholdfilter->Update();

    try 
      {
      typedef itk::ImageFileWriter<MaskImageType> MaskImageFileWriterType;
      MaskImageFileWriterType::Pointer maskwriter = MaskImageFileWriterType::New();
      maskwriter->SetInput(thresholdfilter->GetOutput());
      maskwriter->SetFileName(vm["threshold-mask"].as<std::string>().c_str());
      maskwriter->Update();
      }
    catch (itk::ExceptionObject & e)
      {
      std::cerr << "Could not write threshold mask file" << std::endl;
      std::cerr << e << std::endl;
      }

    }
  
  // Estimate tensors
  typedef itk::ImageToImageFilter<VectorImageType, TensorImageType> DiffusionEstimationBaseType;
  TensorImageType::Pointer tensors;

  if(VERBOSE)
    {
    std::cout << "Estimation method: " << vm["method"].as<EstimationType>() << std::endl;
    }

  DiffusionEstimationFilterType::Pointer llsestimator = DiffusionEstimationFilterType::New();
  
  llsestimator->SetGradientImage(gradientContainer,dwi);
  llsestimator->SetBValue(b0);
  llsestimator->SetThreshold(threshold);
  llsestimator->SetNumberOfThreads(1);
  llsestimator->Update();
  tensors = llsestimator->GetOutput();

  if(vm["method"].as<EstimationType>() == Linear)
    {
    }
  else if(vm["method"].as<EstimationType>() == Nonlinear)
    {
    NLDiffusionEstimationFilterType::Pointer estimator = NLDiffusionEstimationFilterType::New();

    TensorImageType::Pointer llstensors = tensors;

    estimator->SetGradientImage(gradientContainer,dwi);
    estimator->SetBValue(b0);
    estimator->SetThreshold(threshold);
    estimator->SetInitialTensor(llstensors);
    estimator->SetStep(step);
    estimator->SetNumberOfThreads(1);
    estimator->Update();
    tensors = estimator->GetOutput();
    }
  else if(vm["method"].as<EstimationType>() == Weighted)
    {
    WLDiffusionEstimationFilterType::Pointer estimator = WLDiffusionEstimationFilterType::New();

    TensorImageType::Pointer llstensors = tensors;

    estimator->SetGradientImage(gradientContainer,dwi);
    estimator->SetBValue(b0);
    estimator->SetThreshold(threshold);
    estimator->SetInitialTensor(llstensors);
    estimator->SetNumberOfIterations(3);
    estimator->SetNumberOfThreads(1);
    estimator->Update();
    tensors = estimator->GetOutput();
    }
  else if(vm["method"].as<EstimationType>() == MaximumLikelihood)
    {
    MLDiffusionEstimationFilterType::Pointer estimator = MLDiffusionEstimationFilterType::New();
    TensorImageType::Pointer llstensors = tensors;

    estimator->SetGradientImage(gradientContainer,dwi);
    estimator->SetBValue(b0);
    estimator->SetThreshold(threshold);
    estimator->SetInitialTensor(llstensors);
    estimator->SetStep(step);
    estimator->SetNumberOfThreads(1);
    std::cout << "Start sigma: " << sigma << std::endl;
    estimator->SetSigma(sigma);
    estimator->Update();
    tensors = estimator->GetOutput();
    }
  else
    {
    std::cerr << "Invalid estimation method\n"  << std::endl;
    return EXIT_FAILURE;
    }
      
  // wp = D*x
  // wv = M*x
  // wv' = D'M*x

  // Apply measurement frame if it exists
  if(dict.HasKey(NRRD_MEASUREMENT_KEY))
    {
    // measurement frame
    vnl_matrix<double> mf(3,3);
    // imaging frame
    vnl_matrix<double> imgf(3,3);
    
    std::vector<std::vector<double> > nrrdmf;
    itk::ExposeMetaData<std::vector<std::vector<double> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);

    imgf = tensors->GetDirection().GetVnlMatrix();
    for(unsigned int i = 0; i < 3; ++i)
      {
      for(unsigned int j = 0; j < 3; ++j)
        {
        mf(i,j) = nrrdmf[i][j];

        if(i == j)
          nrrdmf[i][j] = 1.0;
        else
          nrrdmf[i][j] = 0.0;
        }
      }

    itk::EncapsulateMetaData<std::vector<std::vector<double> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);

    // prevent slicer error
    if(dict.HasKey("modality"))
      {
      itk::EncapsulateMetaData<std::string>(dict,"modality","DTMRI");
      }


    typedef itk::TensorRotateImageFilter<TensorImageType, TensorImageType, double> TensorRotateFilterType;
    TensorRotateFilterType::Pointer trotate = TensorRotateFilterType::New();
    trotate->SetInput(tensors);
    trotate->SetRotation(vnl_svd<double>(imgf).inverse()*mf);
    trotate->Update();
    tensors = trotate->GetOutput();
    
    }

  // Write tensor file if requested
  try
    {
    typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;

    TensorFileWriterType::Pointer tensorWriter = TensorFileWriterType::New();
    tensorWriter->SetFileName(vm["tensor-output"].as<std::string>().c_str());
    tensors->SetMetaDataDictionary(dict);
    tensorWriter->SetInput(tensors);
    tensorWriter->SetUseCompression(true);
    tensorWriter->Update();
       
    } 
  catch (itk::ExceptionObject e) 
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }   
       
  return EXIT_SUCCESS;
}
