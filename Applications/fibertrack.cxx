/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.4 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "fiberio.h"
#include "pomacros.h"
#include "itkImageToDTIStreamlineTractographyFilter.h"

#include <itkDiffusionTensor3D.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectWriter.h>
#include <itkMetaDataObject.h>

#include <iostream>
#include <string>
#include <cmath>

#include "fibertrackCLP.h"

enum IntegrationType {Euler, Midpoint, RK4};
#if 0
void validate(boost::any& v,
              const std::vector<std::string>& values,
              IntegrationType* target_type,
              int)
{
  using namespace boost::program_options;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s == "euler")
  {
    v = any(Euler);
  }
  else if(s == "modpoint")
  {
    v = any(Midpoint);
  }
  else if(s == "rk4")
  {
    v = any(RK4);
  }
  else
  {
    throw validation_error("Estimation type invalid.  Only \"lls\", \"nls\", \"wls\", and \"ml\" allowed.");
  }
}
#endif

int main(int argc, char* argv[])
{
  typedef itk::DiffusionTensor3D<double> DiffusionTensor;
  typedef itk::Image<DiffusionTensor, 3> TensorImage;
  typedef itk::Image<unsigned short, 3>  LabelImage;
  typedef itk::GroupSpatialObject<3>     FiberBundle;

  typedef itk::ImageFileReader<TensorImage> TensorImageReader;
  typedef itk::ImageFileReader<LabelImage>  LabelImageReader;

  typedef itk::ImageToDTIStreamlineTractographyFilter<TensorImage, LabelImage, FiberBundle> TractographyFilter;

#if 0
  namespace po = boost::program_options;

  po::options_description config("Usage: fibertrack [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("input-tensor-file,i", po::value<std::string>(), "Tensor image")
    ("input-roi-file,r", po::value<std::string>(), "ROI Image")
    ("output-fiber-file,o", po::value<std::string>(), "Fiber file")
    ("source-label,s", po::value<unsigned int>()->default_value(2), "Source label")
    ("target-label,t", po::value<unsigned int>()->default_value(1), "Target label")
    ("forbidden-label,f", po::value<unsigned int>()->default_value(0), "Forbidden label")
    ("max-angle", po::value<double>()->default_value(M_PI/4, "PI/4"), "Maximum angle of change in radians")
    ("step-size", po::value<double>()->default_value(0.5), "Step size in mm")
    ("min-fa", po::value<double>()->default_value(0.2), "Minimum anisotropy")
//    ("integration-method", po::value<RK4>(), "Integration method (euler, midpoint, rk4)")
    ("whole-brain", "Use every voxel in the brain as a potential seed point")
    ("verbose,v", "Verbose output")
    ("really-verbose", "Follow detail of fiber tracking algorithm")
    ("force",
     "Ignore image sanity checks.")
    ;

  po::variables_map vm;
  try
  {
    po::store(po::command_line_parser(argc, argv).
              options(config).run(), vm);
    po::notify(vm);     
  } 
  catch (const po::error &e)
  {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  if(vm.count("help") || !vm.count("input-tensor-file") ||
     !vm.count("input-roi-file") || !vm.count("output-fiber-file"))
  {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
    {
      std::cerr << "Tensor image and roi image needs to be specified." << std::endl;
      return EXIT_FAILURE;
    }
  }
#endif
  PARSE_ARGS;
  
  if(inputTensor == "" || inputROI == "" || outputFiberFile == "")
    {
    std::cerr << "Tensor image and roi image needs to be specified." << std::endl;
    return EXIT_FAILURE;
    }
  TensorImageReader::Pointer tensorreader = TensorImageReader::New();
  LabelImageReader::Pointer  labelreader  = LabelImageReader::New();
  
  tensorreader->SetFileName(inputTensor);
  labelreader->SetFileName(inputROI);

  try
  {
    tensorreader->Update();
    labelreader->Update();
  }
  catch( itk::ExceptionObject & e)
  {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }
  
  if(verbose)
  {
    tensorreader->GetOutput()->Print(std::cout);
  }

  // Sanity check the ROI and tensor image as they must be consistent
  // for the filter to work correctly
  requireequal((tensorreader->GetOutput()->GetSpacing() == labelreader->GetOutput()->GetSpacing()),
               "Image Spacings", force);
  requireequal((tensorreader->GetOutput()->GetLargestPossibleRegion() == labelreader->GetOutput()->GetLargestPossibleRegion()),
               "Image Sizes", force);
  requireequal((tensorreader->GetOutput()->GetOrigin() == labelreader->GetOutput()->GetOrigin()),
               "Image Origins", force);
  requireequal((tensorreader->GetOutput()->GetDirection() == labelreader->GetOutput()->GetDirection()),
               "Image Orientations", force);

  TractographyFilter::Pointer fibertracker = TractographyFilter::New();
  if(reallyVerbose)
    fibertracker->DebugOn();
  if(wholeBrain)
    fibertracker->WholeBrainOn();
  fibertracker->SetTensorImage(tensorreader->GetOutput());
  fibertracker->SetROIImage(labelreader->GetOutput());
  fibertracker->SetSourceLabel(sourceLabel);
  fibertracker->SetTargetLabel(targetLabel);
  fibertracker->SetForbiddenLabel(forbiddenLabel);
  fibertracker->SetMaximumAngleChange(maxAngle);
  fibertracker->SetMinimumFractionalAnisotropy(minFa);
  fibertracker->SetStepSize(stepSize);
  fibertracker->Update();
  
  try
  {
    writeFiberFile(outputFiberFile, fibertracker->GetOutput());
  }
  catch(itk::ExceptionObject e)
  {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
