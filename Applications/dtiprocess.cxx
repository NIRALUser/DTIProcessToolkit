/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2007-11-30 18:44:14 $
  Version:   $Revision: 1.3 $
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

// IO
#include <itkImageFileReader.h>

// dtiprocess headers
#include "transforms.h"
#include "tensoroperations.h"
#include "imageio.h"

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

    // Derived outputs
    ("fa-output,f", po::value<std::string>(), "FA output file")
    ("fa-gradient-output", po::value<std::string>(), "FA gradient")
    ("sigma,s", po::value<double>()->default_value(2.0), "Scale of gradients")
    ("fa-gradmag-output", po::value<std::string>(), "FA gradient magnitude")
    ("color-fa-output,c", po::value<std::string>(), "Color FA output file")
    ("principal-eigenvector-output,V", po::value<std::string>(), "Principal Eigenvector of tensor field")
    ("closest-dotproduct-output,D", po::value<std::string>(),  "Closes dot product of principal eigenvector to all gradient directions")
    ("md-output,m", po::value<std::string>(), "MD output file")
    ("negative-eig-output,n", po::value<std::string>(), "Negative eigenvalue binary image output file")

    // derived output options
    ("scalar-float", "Write scalar [FA,MD] as unscaled float.  Also causes FA to be unscaled [0..1].")
    // ("double", "Writes output in double precisision") // *Currently Default *

    // tensor transformations
    // affine
    ("rot-output,r", po::value<std::string>(),"Rotated tensor output file.  Must also specify the dof file")
    ("dof-file,d", po::value<std::string>(), "DOF for rotation")

    // deformation
    ("deformation-output,w", po::value<std::string>(), "Warped tensor field based on a deformation field.  Must input h as \"h-field\" of transform")
    ("h-field,H", po::value<std::string>(), "HField for warp")
    ("inv-h-field,I", po::value<std::string>(), "Inverse HField for warp")

    // transformation options
    ("interpolation,i", po::value<InterpolationType>()->default_value(Linear,"linear"), "Interpolation type (nearestneighbor, linear, cubic)")
    ("reorientation", po::value<TensorReorientationType>()->default_value(FiniteStrain,"fs (Finite Strain)"), "Reorientation type (fs, ppd)")
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
      return EXIT_SUCCESS;
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

  bool readb0 = false;
  double b0 = 0;

  // Parse gradient directions from image header as specified by the
  // namic conventions defined at http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
  GradientListType::Pointer gradientContainer = GradientListType::New();

  itk::MetaDataDictionary & dict = dtireader->GetOutput()->GetMetaDataDictionary();
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

      if(g.two_norm() > 0.0)
        g = g / g.two_norm();

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

  // Debugging
  if(VERBOSE)
    {
    std::cout << "Interpolation type: " << vm["interpolation"].as<InterpolationType>() << std::endl;
    std::cout << "reorientation type: " << vm["reorientation"].as<TensorReorientationType>() << std::endl;
    }


  double sigma = vm["sigma"].as<double>();

  // Compute FA image
  if(vm.count("fa-output"))
    {
    if(scale)
      writeImage(vm["fa-output"].as<std::string>(), 
                 createFA<unsigned short>(dtireader->GetOutput()));
    else
      writeImage(vm["fa-output"].as<std::string>(), 
                 createFA<double>(dtireader->GetOutput()));
    }

  if(vm.count("fa-gradient-output"))
    {
    writeImage(vm["fa-gradient-output"].as<std::string>(), 
               createFAGradient(dtireader->GetOutput(), sigma));
    }

  if(vm.count("fa-gradmag-output"))
    {
    writeImage(vm["fa-gradmag-output"].as<std::string>(), 
               createFAGradMag(dtireader->GetOutput(), sigma));
    }

  if(vm.count("color-fa-output"))
    {
    writeImage(vm["color-fa-output"].as<std::string>(), 
               createColorFA(dtireader->GetOutput()));
    }

  if(vm.count("principal-eigenvector-output") || 
     vm.count("closest-dotproduct-output"))
    {
    const std::string peo(vm.count("principal-eigenvector-output") ?
                          vm["principal-eigenvector-output"].as<std::string>() :
                          "");

    const std::string cdo(vm.count("closest-dotproduct-output") ?
                          vm["closest-dotproduct-output"].as<std::string>() :
                          "");

    createPrincipalEigenvector(dtireader->GetOutput(),peo,cdo,gradientContainer);
    }

  if(vm.count("md-output"))
    {
    if(scale)
      writeImage(vm["md-output"].as<std::string>(), 
                 createMD<unsigned short>(dtireader->GetOutput()));
    else
      writeImage(vm["md-output"].as<std::string>(), 
                 createMD<double>(dtireader->GetOutput()));
    }

  if(vm.count("negative-eig-output"))
    {
    writeImage(vm["negative-eig-output"].as<std::string>(), 
               createNegativeEigenValueLabel(dtireader->GetOutput()));
    }


  if(vm.count("rot-output"))
    {
    if(!vm.count("dof-file"))
      {
      std::cerr << "Tensor rotation requested, but dof file not specified" << std::endl;
      return EXIT_FAILURE;
      }
    writeImage(vm["rot-output"].as<std::string>(),
               createROT(dtireader->GetOutput(),
                         vm["dof-file"].as<std::string>()));
    }
   

  if(vm.count("deformation-output"))
    {
    if(!vm.count("h-field") || !vm.count("inv-h-field"))
      {
      std::cerr << "Deformation field info not fully specified" << std::endl;
      return EXIT_FAILURE;
      }
    writeImage(vm["deformation-output"].as<std::string>(),
               createWarp(dtireader->GetOutput(),
                          vm["h-field"].as<std::string>(),
                          vm["inv-h-field"].as<std::string>(),
                          vm["reorientation"].as<TensorReorientationType>(),
                          vm["interpolation"].as<InterpolationType>()));
    }

  return EXIT_SUCCESS;
}


