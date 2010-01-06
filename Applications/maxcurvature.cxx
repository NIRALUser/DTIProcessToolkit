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
// STL includes
#include <string>
#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkNthElementImageAdaptor.h>
#include <itkAbsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkVersion.h>
#include "maxcurvatureCLP.h"

enum CurvatureType { MaxEigenvalue, SmoothNormalized, RawNormalized, UnPossible };
#if 0
void validate(boost::any& v,
              const std::vector<std::string>& values,
              CurvatureType* curvature_type,
              int)
{
  using namespace po;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s ==  "orig")
    {
    v = any(MaxEigenvalue);
    }
  else if (s == "snorm")
    {
    v = any(SmoothNormalized);
    } 
  else if (s == "rnorm")
    {
    v = any(RawNormalized);
    }
  else
    {
    throw validation_error("Curvature type invalid.  Only \"orig\", \"snorm\", and \"rnorm\" allowed.");
    }


}
#endif

int main(int argc, char* argv[])
{
#if 0
  // Read program options/configuration
  po::options_description config("Usage: maxcurvature input-image [options]");
  config.add_options()
    // General options
    ("help,h", "produce this help message")
    ("verbose,v","produces verbose output")

    // Outputs
    ("output,o", po::value<std::string>(), "Output file")
    ("sigma,s", po::value<double>()->default_value(2.0), "Scale of gradients")
    ("type,t", po::value<CurvatureType>()->default_value(MaxEigenvalue,"orig"), "Curvature type")
    ;
  
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("image", po::value<std::string>(), "FA image")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("image",1);


  po::variables_map vm;

  try
    {
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);     
    } 
  catch (const po::error &e)
    {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
    }

  // End option reading configuration
#endif
  PARSE_ARGS;
  
  // Display help if asked or program improperly called
  if(image == "")
    {
    std::cerr << "maxcurvature: Missing input FA image (--image)" << std::endl;
    return EXIT_FAILURE;
    }

  // bool VERBOSE = verbose;

  if(output == "")
    {
    std::cerr << "maxcurvature: Must specify output file" << std::endl;
    return EXIT_FAILURE;
    }

  CurvatureType ctype =
    type == "orig" ? MaxEigenvalue :
    (type == "snorm" ? SmoothNormalized :
     (type == "rnorm" ? RawNormalized : UnPossible));
  if(ctype == UnPossible)
    {
    std::cerr << "maxcurvature: Bad Curvature type " << type << std::endl;
    return EXIT_FAILURE;
    }

  typedef unsigned short PixelType;
  typedef double FloatPixelType;
  const int DIM = 3;
  typedef itk::SymmetricSecondRankTensor<double, DIM> HessianPixelType;
  typedef itk::Vector<double, DIM> EigenValuePixelType;

  typedef itk::Image<PixelType, DIM> ImageType;
  typedef itk::Image<FloatPixelType, DIM> FloatImageType;

  typedef itk::Image<HessianPixelType, DIM> HessianImageType;
  typedef itk::Image<EigenValuePixelType, DIM> EigenValueImageType;

  typedef itk::ImageFileReader<ImageType> FileReaderType;
 
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(image);

  // sigma set by PARSE_ARGS
  //  double sigma = vm["sigma"].as<double>();
  typedef itk::HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;
  HessianFilterType::Pointer hessian = HessianFilterType::New();
  hessian->SetInput(reader->GetOutput());
  hessian->SetSigma(sigma);

  typedef itk::SymmetricEigenAnalysisImageFilter<HessianImageType,EigenValueImageType> EigenAnalysisFilterType;
  EigenAnalysisFilterType::Pointer eigen = EigenAnalysisFilterType::New();
  eigen->SetInput(hessian->GetOutput());
  eigen->OrderEigenValuesBy(EigenAnalysisFilterType::FunctorType::OrderByValue);
  eigen->SetDimension(3);

  eigen->Update();
  typedef itk::NthElementImageAdaptor<EigenValueImageType,FloatPixelType> ElementSelectAdaptorType;
  ElementSelectAdaptorType::Pointer elementSelect = ElementSelectAdaptorType::New();
  elementSelect->SetImage(eigen->GetOutput());
  elementSelect->SelectNthElement(0);

  typedef itk::CastImageFilter<ElementSelectAdaptorType,FloatImageType> CastFilter2Type;
  CastFilter2Type::Pointer cast2 = CastFilter2Type::New();
  cast2->SetInput(elementSelect);
  cast2->Update();
  
  typedef itk::ShiftScaleImageFilter<FloatImageType,FloatImageType> ScaleImageType;
  ScaleImageType::Pointer scale = ScaleImageType::New();
  scale->SetShift(0.0);

  if(ctype == SmoothNormalized)
    {
    typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, FloatImageType> SmoothType;
    SmoothType::Pointer smooth = SmoothType::New();
    smooth->SetInput(reader->GetOutput());
    smooth->SetSigma(sigma);
    
    typedef itk::DivideImageFilter<FloatImageType, FloatImageType, FloatImageType> DivideType;
    DivideType::Pointer divider = DivideType::New();
    divider->SetInput1(cast2->GetOutput());
    divider->SetInput2(smooth->GetOutput());

    scale->SetScale(-10000.0);
    scale->SetInput(divider->GetOutput());
    }
  else if(ctype == MaxEigenvalue)
    {
    scale->SetScale(-10.0);
    scale->SetInput(cast2->GetOutput());
    }

  typedef itk::IntensityWindowingImageFilter<FloatImageType> WindowFilterType;
  WindowFilterType::Pointer window = WindowFilterType::New();
  window->SetInput(scale->GetOutput());
  window->SetWindowMinimum(0);
  window->SetWindowMaximum(32768);
  window->SetOutputMinimum(0);
  window->SetOutputMaximum(32768);

  typedef itk::CastImageFilter<FloatImageType,ImageType> CastFilterType;
  CastFilterType::Pointer cast = CastFilterType::New();
  cast->SetInput(window->GetOutput());

  typedef itk::ImageFileWriter<ImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetUseCompression(true);
  writer->SetInput(cast->GetOutput());
  writer->SetFileName(output);

  try
    {
    eigen->Update();
    writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e << std::endl;
    return -1;
    }


  return 0;
}
