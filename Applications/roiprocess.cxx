/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $
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
#include <fstream>

// boost includes
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/cmdline.hpp>

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkWarpImageFilter.h>

#include <itkImageDuplicator.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkAffineTransform.h>

enum InterpolationType {NearestNeighbor, Linear};

// Validates the interpolation type option string to the the allowed
// values for interpolation methods.  Currently nearestneighbor,
// linear, or cubic.
void validate(boost::any& v,
              const std::vector<std::string>& values,
              InterpolationType* target_type,
              int)
{
  using namespace boost::program_options;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s ==  "nearestneighbor")
    {
    v = any(NearestNeighbor);
    }
  else if (s == "linear")
    {
    v = any(Linear);
    } 
  else
    {
    throw validation_error("Interpolation type invalid.  Only \"nearestneighbor\" and \"linear\"\"cubic\" allowed.");
    }
}

#include "itkHFieldToDeformationFieldImageFilter.h"
int main(int argc, char* argv[])
{
  namespace po = boost::program_options;

  // Read program options/configuration
  po::options_description config("Usage: roiprocess input-roi h-field output-roi [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v","produces verbose output")
    ("interpolation,i", po::value<InterpolationType>()->default_value(NearestNeighbor, "nearestneighbor"), "Interpolation type (nearestneighbor, linear)")
    ("threshold,t", po::value<double>()->default_value(0.1), "Threshold for linear interpolation")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("roi-file", po::value<std::string>(), "Label file")
    ("h-field", po::value<std::string>(), "HField for warp")
    ("roi-output", po::value<std::string>(), "Warped roi file based on a deformation field.")
    ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("roi-file",1);
  p.add("h-field",1);
  p.add("roi-output",1);

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
  
  if(vm.count("help") || !vm.count("roi-file") ||
     !vm.count("roi-output") || !vm.count("h-field"))
  {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      return EXIT_FAILURE;
  }

  // Reading roi image
  std::string roifile(vm["roi-file"].as<std::string>());
  std::string warpfile(vm["h-field"].as<std::string>());
  std::string outputfile(vm["roi-output"].as<std::string>());

  typedef unsigned char ROIPixelType;
  typedef itk::Image<ROIPixelType, 3> ROIImageType;
  typedef itk::ImageFileReader<ROIImageType> ROIImageReaderType;
  ROIImageReaderType::Pointer roireader = ROIImageReaderType::New();
  roireader->SetFileName(roifile);
  roireader->Update();

  typedef itk::ImageDuplicator<ROIImageType> DuplicatorType;
  DuplicatorType::Pointer duplicate = DuplicatorType::New();
  duplicate->SetInputImage(roireader->GetOutput());
  duplicate->Update();

  // Reading deformation field
  typedef itk::Vector<double, 3> DeformationPixelType;
  typedef itk::Image<DeformationPixelType, 3> DeformationImageType;
  typedef itk::ImageFileReader<DeformationImageType> DeformationImageReader;
  DeformationImageReader::Pointer defreader = DeformationImageReader::New();
  defreader->SetFileName(warpfile);
  defreader->Update();
  
  typedef itk::LinearInterpolateImageFunction<ROIImageType, double> DeformationInterpolateType;
  DeformationInterpolateType::Pointer definterp = DeformationInterpolateType::New();
  definterp->SetInputImage(roireader->GetOutput());

  typedef itk::NearestNeighborInterpolateImageFunction<ROIImageType, double> ROIInterpolateType;
  ROIInterpolateType::Pointer nninterp = ROIInterpolateType::New();
  nninterp->SetInputImage(roireader->GetOutput());

  ROIImageType::Pointer output = duplicate->GetOutput();
  // iterator over duplicated point
  typedef itk::ImageRegionIterator<ROIImageType> ROIIteratorType;
  typedef itk::ImageRegionConstIterator<DeformationImageType> HFieldIteratorType;
  ROIIteratorType outputIt(output,output->GetLargestPossibleRegion());
  HFieldIteratorType hfieldIt(defreader->GetOutput(),
                              defreader->GetOutput()->GetLargestPossibleRegion());

  outputIt.GoToBegin();
  hfieldIt.GoToBegin();
  for( ; !outputIt.IsAtEnd(); ++outputIt, ++hfieldIt)
  {
    ROIInterpolateType::ContinuousIndexType ci;
    DeformationImageType::PixelType h = hfieldIt.Get();
    ci[0] = h[0];
    ci[1] = h[1];
    ci[2] = h[2];
    
    ROIPixelType result = 0;
    if(!output->GetLargestPossibleRegion().IsInside(ci))
      result = 0;
    else if(vm["interpolation"].as<InterpolationType>() == NearestNeighbor)
      result = static_cast<ROIPixelType>(nninterp->EvaluateAtContinuousIndex(ci));
    else if(vm["interpolation"].as<InterpolationType>() == Linear)
      result = definterp->EvaluateAtContinuousIndex(ci) > vm["threshold"].as<double>() ? 1 : 0;
    else
      throw std::logic_error("Invalid interpolation type");
  }
  if(!hfieldIt.IsAtEnd())
  {
    std::cerr << "Image and hfield of different size" << std::endl;
  }

  typedef itk::ImageFileWriter<ROIImageType> ImageWriterType;
  ImageWriterType::Pointer writer =  ImageWriterType::New();
  writer->SetFileName(outputfile);
  writer->SetInput(output);
  writer->Update();

  return EXIT_SUCCESS;
}
