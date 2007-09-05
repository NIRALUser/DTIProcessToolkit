/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2007-09-05 19:35:36 $
  Version:   $Revision: 1.2 $
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

#include "itkHFieldToDeformationFieldImageFilter.h"
int main(int argc, char* argv[])
{
  namespace po = boost::program_options;

  // Read program options/configuration
  po::options_description config("Usage: roiprocess input-roi output-roi [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v","produces verbose output")
    ("h-field,H", po::value<std::string>(), "HField for warp")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("roi-file", po::value<std::string>(), "DTI roi file")
    ("roi-output", po::value<std::string>(), "Warped roi file based on a deformation field.  Must input h as \"h-field\" of transform")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("roi-file",1);
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
  roireader->SetFileName(roifile.c_str());
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
  defreader->SetFileName(warpfile.c_str());
  defreader->Update();

//   typedef itk::LinearInterpolateImageFunction<DeformationImageType,double> DeformationInterpolateType;
//   DeformationInterpolateType::Pointer definterp = DeformationInterpolateType::New();
//   definterp->SetInputImage(defconv->GetOutput());

  typedef itk::NearestNeighborInterpolateImageFunction<ROIImageType,double> ROIInterpolateType;
  ROIInterpolateType::Pointer nninterp = ROIInterpolateType::New();
  nninterp->SetInputImage(roireader->GetOutput());
//  nninterp->

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
    
    if(!output->GetLargestPossibleRegion().IsInside(ci))
      outputIt.Set(0);
    else
      outputIt.Set(static_cast<unsigned char>(nninterp->EvaluateAtContinuousIndex(ci)));
    }
  if(!hfieldIt.IsAtEnd())
    {
    std::cerr << "Image and hfield of different size" << std::endl;
    }

  typedef itk::ImageFileWriter<ROIImageType> ImageWriterType;
  ImageWriterType::Pointer writer =  ImageWriterType::New();
  writer->SetFileName(outputfile.c_str());
  writer->SetInput(output);
  writer->Update();

  return EXIT_SUCCESS;
}
