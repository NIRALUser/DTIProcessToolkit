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
#include "itkVectorMaskNegatedImageFilter.h"
#include "itkVectorMaskImageFilter.h"
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>

#include <itkNthElementImageAdaptor.h>

namespace po = boost::program_options;

#include <itkTensorModelResidualImageFilter.h>

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
  
  //  unsigned int scale;
  
  // Read program options/configuration
  po::options_description config("Usage: dtiestim dwi-image residual-output [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("brain-mask,M", po::value<std::string>(), "sets brain mask")
    ("bad-region-mask,B", po::value<std::string>(), "sets bad region mask")
    ("double,d", "Writes tensors in double precision -- not implemented yet")
    ("threshold,t", po::value<DWIPixelType>(),"Baseline threshold for estimation")
    ("sigma,s", po::value<double>(), "Rician noise parameter")
    ("verbose,v", "Verbose output")
  ;
//    ("fa-scale,s", po::value<unsigned int>(&scale)->default_value(10000),"FA scale factor.  If set the FA value is scaled by this factor and written out in an integer image format.")

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("dwi-image", po::value<std::string>(), "DWI image volume")
    ("residual-output", po::value<std::string>(), "Residual output")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("dwi-image",1);
  p.add("residual-output",1);

  po::variables_map vm;

  try
    {
    po::store(po::command_line_parser(argc, argv).
            options(all).positional(p).run(), vm);
    po::notify(vm);     
    } 
  catch (const po::error &e)
    {
    //std::cerr << e << std::endl;
    std::cout << config << std::endl;
    return EXIT_FAILURE;
    }

  // End option reading configuration

  // Display help if asked or program improperly called
  if(vm.count("help") || !vm.count("dwi-image") || !vm.count("residual-output"))
    {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      return EXIT_FAILURE;
    }

  bool VERBOSE = false;
  if(vm.count("verbose"))
    VERBOSE = true;

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
  typedef itk::DiffusionTensor3DReconstructionImageFilter<DWIPixelType,DWIPixelType, RealType>
    DiffusionEstimationFilterType;
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

  if(vm.count("brain-mask") && vm.count("bad-region-mask"))
    {
    std::cerr << "--brain-mask and --bad-region-mask are mutually exclusive" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileReader<MaskImageType> MaskFileReaderType;
  MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();

  // Read brain mask
  if(vm.count("brain-mask"))
    {

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
  
  // Read negatice mask
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

  // If we didnt specify a threshold compute it as the mean of the
  // baseline image
  if(!vm.count("threshold"))
    {
    

    return EXIT_FAILURE;
    }
  
  // Estimate tensors
  DiffusionEstimationFilterType::Pointer estimator = DiffusionEstimationFilterType::New();
  estimator->SetGradientImage(gradientContainer,dwi);
  estimator->SetBValue(b0);
  estimator->SetNumberOfThreads(1);
  estimator->SetThreshold(vm["threshold"].as<DWIPixelType>());
  estimator->Update();
     
  typedef DiffusionEstimationFilterType::OutputImageType TensorImageType;

  typedef itk::TensorModelResidualImageFilter<VectorImageType,
    TensorImageType, RealImageType> ModelResidualFilterType;

  ModelResidualFilterType::Pointer residual = ModelResidualFilterType::New();
  residual->SetInput1(dwi);
  residual->SetInput2(estimator->GetOutput());
  residual->SetBValue(b0);
  residual->SetSigma(vm["sigma"].as<double>());
  residual->SetGradientList(gradientContainer);
  residual->Update();

  // Write tensor file if requested
  try
    {
    typedef itk::ImageFileWriter<RealImageType> ResidualFileWriterType;

    ResidualFileWriterType::Pointer residualWriter = ResidualFileWriterType::New();
    residualWriter->SetFileName(vm["residual-output"].as<std::string>().c_str());
    residualWriter->SetInput(residual->GetOutput());
    residualWriter->Update();
       
    } 
  catch (itk::ExceptionObject e) 
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }   
       
  return EXIT_SUCCESS;
}
