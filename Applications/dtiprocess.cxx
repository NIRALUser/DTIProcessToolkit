/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009-12-14 12:29:41 $
  Version:   $Revision: 1.8 $
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

// boost includes
#include <boost/program_options.hpp>

// ITK includes
// datastructures
#include <itkMetaDataObject.h>
#include <itkVersion.h>

// IO
#include <itkImageFileReader.h>

// dtiprocess headers
#include "itkVectorMaskImageFilter.h"
#include "transforms.h"
#include "tensoroperations.h"
#include "imageio.h"
#include "deformationfieldio.h"

namespace po = boost::program_options;

// Bad global variables.  TODO: remove these
bool VERBOSE=false;

// Validates the interpolation type option string to the the allowed
// values for interpolation methods.  Currently nearestneighbor,
// linear, or cubic.
void validate(boost::any& v,
              const std::vector<std::string>& values,
              InterpolationType* target_type,
              int)
{
  using namespace po;
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
  else if (s == "cubic")
  {
    v = any(Cubic);
  }
  else
  {
    throw validation_error("Interpolation type invalid.  Only \"nearestneighbor\", \"linear\", and \"cubic\" allowed.");
  }
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              TensorReorientationType* target_type,
              int)
{
  using namespace po;
  using boost::any;
  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s == "fs")
  {
    v = any(FiniteStrain);
  }
  else if(s == "ppd")
  {
    v = any(PreservationPrincipalDirection);
  }
  else
  {
    throw validation_error("Reorientation type invalid.  Only \"fs\" or \"ppd\"");
  }
}


int main(int argc, char* argv[])
{

  // Read program options/configuration
  po::options_description config("Usage: dtiprocess tensor-image [options]");
  config.add_options()
    // General options
    ("help,h", "produce this help message")
    ("verbose,v","produces verbose output")
    ("mask", po::value<std::string>(), "Masks tensors")

    // Derived outputs
    ("fa-output,f", po::value<std::string>(), "FA output file")
    ("md-output,m", po::value<std::string>(), "MD output file")
    ("fa-gradient-output", po::value<std::string>(), "FA gradient")
    ("sigma,s", po::value<double>()->default_value(2.0), "Scale of gradients.")
    ("fa-gradmag-output", po::value<std::string>(), "FA gradient magnitude")
    ("color-fa-output,c", po::value<std::string>(), "Color FA output file")
    ("principal-eigenvector-output,V", po::value<std::string>(), "Principal Eigenvector of tensor field")
//    ("closest-dotproduct-output,D", po::value<std::string>(),  "Closes dot product of principal eigenvector to all gradient directions")
    ("negative-eig-output,n", po::value<std::string>(), "Negative eigenvalue binary image output file")
    ("frobenius-norm-output", po::value<std::string>(), "Frobenius norm output")

    ("lambda1-output", po::value<std::string>(), "Lambda 1 (largest eigenvalue) output")
    ("lambda2-output", po::value<std::string>(), "Lambda 2 (middle eigenvalue) output")
    ("lambda3-output", po::value<std::string>(), "Lambda 3 (smallest eigenvalue) output")
  

    // derived output options
    ("scalar-float", "Write scalar [FA,MD] as unscaled float.  Also causes FA to be unscaled [0..1].")
    // ("double", "Writes output in double precisision") // *Currently Default *

    // tensor transformations
    // affine
    ("rot-output,r", po::value<std::string>(),"Rotated tensor output file.  Must also specify the dof file.")
    ("dof-file,d", po::value<std::string>(),"Transformation file for affine transformation. This must be RView text format (not binary!).") //This can be RView or ITK format.")

    // deformation
    ("deformation-output,w", po::value<std::string>(), "Warped tensor field based on a deformation field.  This option requires the --forward,-F transformation to be specified.")
    ("forward,F", po::value<std::string>(), "Forward transformation.  Assumed to be a deformation field in world coordinates, unless the --h-field option is specified.")
    //("inverse,I", po::value<std::string>(), "Inverse of warp (DEPRECATED: NO LONGER REQUIRED)")
    ("h-field", "forward and inverse transformations are h-fields instead of displacement fields")

    // transformation options
    ("interpolation,i", po::value<InterpolationType>()->default_value(Linear,"linear"), "Interpolation type (nearestneighbor, linear, cubic)")
    ("reorientation", po::value<TensorReorientationType>()->default_value(FiniteStrain,"fs (Finite Strain)"), "Reorientation type (fs, ppd)")

    // statistics options
    // ("stats", po::value<std::string>(), "Compute statistics")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("dti-image", po::value<std::string>(), "DTI tensor volume")
    ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("dti-image",1);

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

  // Display help if asked or program improperly called
  if(vm.count("help") || !vm.count("dti-image"))
  {
    std::cout << config << std::endl;
    if(vm.count("help"))
    {
      std::cout << "$Date: 2009-12-14 12:29:41 $ $Revision: 1.8 $" << std::endl;
      std::cout << ITK_SOURCE_VERSION << std::endl;
      return EXIT_SUCCESS;
    }
    else
      return EXIT_FAILURE;
  }

  // If the value scale is true (default) we scal FA and MD values to
  // integer ranges.
  bool scale = true;
  if(vm.count("scalar-float"))
  {
    scale = false;
  }

  if(vm.count("verbose"))
  {
    VERBOSE = true;
  }

  // Read tensor image
  typedef itk::ImageFileReader<TensorImageType> FileReaderType;
  FileReaderType::Pointer dtireader = FileReaderType::New();
  dtireader->SetFileName(vm["dti-image"].as<std::string>().c_str());

  try
  {
    dtireader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << e <<std::endl;
    return EXIT_FAILURE;
  }

  // Parse gradient directions from image header as specified by the
  // namic conventions defined at http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
  GradientListType::Pointer gradientContainer = GradientListType::New();

  itk::MetaDataDictionary & dict = dtireader->GetOutput()->GetMetaDataDictionary();


  std::vector<std::string> keys = dict.GetKeys();
  for(std::vector<std::string>::const_iterator it = keys.begin();
      it != keys.end(); ++it)
  {
    if( it->find("DWMRI_gradient") != std::string::npos)
    {
      std::string value;

      itk::ExposeMetaData<std::string>(dict, *it, value);
      std::istringstream iss(value);
      GradientType g;
      iss >> g[0] >> g[1] >> g[2];

      unsigned int ind;
      std::string temp = it->substr(it->find_last_of('_')+1);
      ind = atoi(temp.c_str());
      
      gradientContainer->InsertElement(ind,g);
    }
  }
  for(std::vector<std::string>::const_iterator it = keys.begin();
      it != keys.end(); ++it)
  {
    if( it->find("DWMRI_NEX") != std::string::npos)
    {
      std::string numrepstr;

      itk::ExposeMetaData<std::string>(dict, *it, numrepstr);
      unsigned int numreps = atoi(numrepstr.c_str());

      std::string indtorepstr = it->substr(it->find_last_of('_')+1);
      unsigned int indtorep =  atoi(indtorepstr.c_str());

      GradientType g = gradientContainer->GetElement(indtorep);

      for(unsigned int i = indtorep+1; i < indtorep+numreps; i++)
        gradientContainer->InsertElement(i,g);
    }

  }

  // Debugging
  if(VERBOSE)
  {
    std::cout << "Interpolation type: " << vm["interpolation"].as<InterpolationType>() << std::endl;
    std::cout << "reorientation type: " << vm["reorientation"].as<TensorReorientationType>() << std::endl;
  }

  TensorImageType::Pointer tensors = dtireader->GetOutput();
  if(vm.count("mask"))
  {
    typedef itk::ImageFileReader<LabelImageType> MaskFileReaderType;
    MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();
    maskreader->SetFileName(vm["mask"].as<std::string>());

    typedef itk::VectorMaskImageFilter<TensorImageType,LabelImageType,TensorImageType> MaskFilterType;
    MaskFilterType::Pointer mask = MaskFilterType::New();
    mask->SetInput1(tensors);
    mask->SetInput2(maskreader->GetOutput());

    try
    {
      mask->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << e <<std::endl;
      return EXIT_FAILURE;
    }

    tensors = mask->GetOutput();
  }

  double sigma = vm["sigma"].as<double>();

  // Compute FA image
  if(vm.count("fa-output"))
  {
    if(scale)
      writeImage(vm["fa-output"].as<std::string>(), 
                 createFA<unsigned short>(tensors));
    else
      writeImage(vm["fa-output"].as<std::string>(), 
                 createFA<double>(tensors));
  }


  if(vm.count("fa-gradient-output"))
  {
    writeImage(vm["fa-gradient-output"].as<std::string>(), 
               createFAGradient(tensors, sigma));
  }

  if(vm.count("fa-gradmag-output"))
  {
    writeImage(vm["fa-gradmag-output"].as<std::string>(), 
               createFAGradMag(tensors, sigma));
  }

  if(vm.count("color-fa-output"))
  {
    writeImage(vm["color-fa-output"].as<std::string>(), 
               createColorFA(tensors));
  }

  if(vm.count("principal-eigenvector-output"))// || 
//     vm.count("closest-dotproduct-output"))
  {
    writeImage(vm["principal-eigenvector-output"].as<std::string>(),
               createPrincipalEigenvector(tensors));
  }

  if(vm.count("md-output"))
  {
    if(scale)
      writeImage(vm["md-output"].as<std::string>(), 
                 createMD<unsigned short>(tensors));
    else
      writeImage(vm["md-output"].as<std::string>(), 
                 createMD<double>(tensors));
  }

  if(vm.count("lambda1-output"))
  {
    if(scale)
      writeImage(vm["lambda1-output"].as<std::string>(), 
                 createLambda<unsigned short>(tensors, Lambda1));
    else
      writeImage(vm["lambda1-output"].as<std::string>(), 
                 createLambda<double>(tensors, Lambda1));
  }

  if(vm.count("lambda2-output"))
  {
    if(scale)
      writeImage(vm["lambda2-output"].as<std::string>(), 
                 createLambda<unsigned short>(tensors, Lambda2));
    else
      writeImage(vm["lambda2-output"].as<std::string>(), 
                 createLambda<double>(tensors, Lambda2));
  }

  if(vm.count("lambda3-output"))
  {
    if(scale)
      writeImage(vm["lambda3-output"].as<std::string>(), 
                 createLambda<unsigned short>(tensors, Lambda3));
    else
      writeImage(vm["lambda3-output"].as<std::string>(), 
                 createLambda<double>(tensors, Lambda3));
  }

  if(vm.count("frobenius-norm-output"))
  {
    if(scale)
      writeImage(vm["frobenius-norm-output"].as<std::string>(),
                 createFro<unsigned short>(tensors));
    else
      writeImage(vm["frobenius-norm-output"].as<std::string>(),
                 createFro<double>(tensors));
  }

  if(vm.count("negative-eig-output"))
  {
    writeImage(vm["negative-eig-output"].as<std::string>(), 
               createNegativeEigenValueLabel(tensors));
  }

  if(vm.count("rot-output"))
  {
    if(!vm.count("dof-file"))
    {
      std::cerr << "Tensor rotation requested, but dof file not specified" << std::endl;
      return EXIT_FAILURE;
    }
    writeImage(vm["rot-output"].as<std::string>(),
               createROT(tensors,
                         vm["dof-file"].as<std::string>()));
  }
   

  if(vm.count("deformation-output"))
  {
    if(!vm.count("forward"))
    {
      std::cerr << "Deformation field info not fully specified" << std::endl;
      return EXIT_FAILURE;
    }
    DeformationImageType::Pointer forward;
    
    DeformationFieldType dftype = Displacement;
    if(vm.count("h-field"))
      dftype = HField;
    
    forward = readDeformationField(vm["forward"].as<std::string>(), dftype);

    writeImage(vm["deformation-output"].as<std::string>(),
               createWarp(tensors,
                          forward,
                          vm["reorientation"].as<TensorReorientationType>(),
                          vm["interpolation"].as<InterpolationType>()));
  }

  if(vm.count("stats"))
  {
    if(!vm.count("mask"))
    {
      std::cerr << "WARNING: No mask specified.  Computing whole brain statistics" << std::endl;
    }
    if(VERBOSE)
      std::cout << "Computing tensor stats" << std::endl;
    
//    TensorRegionStatistics tstat = computeTensorStatistics(tensors);
    
  }

  return EXIT_SUCCESS;
}


