/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>
#include <string>
#include <vector>

#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVector.h>
#include <itkExtractImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkVersion.h>

#include "itkLogEuclideanTensorImageFilter.h"
#include "itkExpEuclideanTensorImageFilter.h"
#include "dtiaverageCLP.h"

template <class TElementType>
class PixelDivider
{
public:
  // Divide by zero so we get warned if we didnt replace the default functor
  PixelDivider() : m_Denominator(0) {}
  PixelDivider(double denominator) : m_Denominator(denominator) {}

  TElementType operator()(const TElementType &numerator) { return numerator / m_Denominator; }

  bool operator !=(const PixelDivider & rhs) { return true; }

  double m_Denominator;
};

enum StatisticsType {Euclidean, LogEuclidean, PGA};
#if 0
void validate(boost::any& v,
              const std::vector<std::string>& values,
              StatisticsType* target_type,
              int)
{
  using namespace boost::program_options;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s == "euclid" || s == "euclidean")
  {
    v = any(Euclidean);
  }
  else if (s == "le" || s == "log-euclidean")
  {
    v = any(LogEuclidean);
  }
  else if (s == "pga")
  {
    v = any(PGA);
  }
  else
  {
    throw validation_error("Statistics type invalid.  Only \"euclidean\", \"log-euclidean\", and \"pga\" allowed.");
  }
}
#endif

int main(int argc, char* argv[])
{
#if 0
  namespace po = boost::program_options;

  // Read program options/configuration
  po::options_description config("Usage: dtiaverage tensor-output input1 [input2] [input3] [...] [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("method,m", po::value<StatisticsType>()->default_value(Euclidean,"Euclidean"),
     "Statistics method (euclidean,log-euclidean,pga)")
    ("verbose,v",
     "Verbose output")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("tensor-output", po::value<std::string>(), "Averaged tensor volume.")
    ("inputs", po::value<std::vector<std::string> >(), "Tensor inputs.")
    ;
  
  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("tensor-output",1);
  p.add("inputs",-1);

  po::variables_map vm;

  try
  {
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);     
  } 
  catch (const po::error &e)
  {
    std::cerr << "Error parsing arguments" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << config << std::endl;
    return EXIT_FAILURE;
  }

  if(vm.count("help"))
  {
    std::cout << config << std::endl;
    std::cout << "Version: $Date: 2009/01/09 15:39:51 $ $Revision: 1.5 $" << std::endl;
    std::cout << ITK_SOURCE_VERSION << std::endl;
    return EXIT_SUCCESS;
  }
#endif
  PARSE_ARGS;
  
  typedef double RealType;
  typedef itk::DiffusionTensor3D<RealType> TensorPixelType;
  typedef itk::Image<TensorPixelType, 3> TensorImageType;
  typedef itk::ImageFileReader<TensorImageType> TensorFileReader;
  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType LogTensorImageType;
  typedef LogTensorImageType::PixelType LogTensorPixelType;

  typedef itk::AddImageFilter<LogTensorImageType, LogTensorImageType, LogTensorImageType> AddImageFilter;
  typedef itk::UnaryFunctorImageFilter<LogTensorImageType, LogTensorImageType, PixelDivider<LogTensorPixelType> > DivideImageFilter;

  typedef itk::ImageDuplicator<LogTensorImageType> DuplicateImageFilter;

  //  const std::vector<std::string> sources = vm["inputs"].as<std::vector<std::string> >();
  
  const int n = inputs.size();
  
  TensorFileReader::Pointer reader = TensorFileReader::New();
  AddImageFilter::Pointer adder = AddImageFilter::New();
  LogEuclideanFilter::Pointer logfilt = LogEuclideanFilter::New();

  DuplicateImageFilter::Pointer dup = DuplicateImageFilter::New();

  if(verbose)
    std::cout << "Loading: " <<  inputs[0] << std::endl;
  reader->SetFileName(inputs[0]);
  reader->Update();
  logfilt->SetInput(reader->GetOutput());
  logfilt->Update();
 
  dup->SetInputImage(logfilt->GetOutput());
  dup->Update();
  LogTensorImageType::Pointer average = dup->GetOutput();

  for(int i = 1; i < n; ++i)
  {
  if(verbose)
      std::cout << "Loading: " <<  inputs[i] << std::endl;
    reader->SetFileName(inputs[i].c_str());
    
    adder->SetInput1(average);
    adder->SetInput2(logfilt->GetOutput());
    adder->Update();
    
    dup->SetInputImage(adder->GetOutput());
    dup->Update();
    average = dup->GetOutput();
    
  }

  DivideImageFilter::Pointer divide = DivideImageFilter::New();
  divide->SetFunctor(PixelDivider<LogTensorPixelType>(n));
  divide->SetInput(average);
  divide->Update();

  typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
  ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
  expf->SetInput(divide->GetOutput());
  
  typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;
  TensorFileWriterType::Pointer twrit = TensorFileWriterType::New();
  twrit->SetUseCompression(true);
  twrit->SetInput(expf->GetOutput());
  twrit->SetFileName(tensorOutput);
  twrit->Update();
  
  return EXIT_SUCCESS;
}

