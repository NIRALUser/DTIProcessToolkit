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
#include <itkDiffusionTensor3D.h>

#include <itkImageFileReader.h>

#include <itkVectorLinearInterpolateImageFunction.h>
#include "itkVectorBSplineInterpolateImageFunction.h"

#include "FiberCalculator.h"
#include "deformationfieldoperations.h"
#include "fiberio.h"
#include "dtitypes.h"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  // Read program options/configuration
  po::options_description config("Usage: fiberprocess input-fiber [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v", "produces verbose output")
    ("fiber-output,o", po::value<std::string>(), "Output fiber file.  May be warped or updated with new data depending on other options used.")

    ("h-field,H", po::value<std::string>(), "HField for warp and statistics lookup.  If this option is used tensor-volume must also be specified.")
    ("no-warp,n", "Do not warp the geometry of the tensors only obtain the new statistics")
    ("tensor-volume,T", po::value<std::string>(), "Interpolate tensor values from the given field")

    // ******  TODO **********
    // Not yet implemented
    // ("voxelize,V", po::value<std::string>(),"Voxelize fiber into a label map.")
    // ("voxel-label,l", po::value<unsigned int>()->default_value(1),"Label for voxelized fiber")

    // Compute 1-D statistics
    // ("mean-statistics,s", po::value<std::string>(), "Write summary statistics to text file")
    // ******  TODO **********
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("fiber-file", po::value<std::string>(), "DTI fiber file")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("fiber-file",1);

  po::variables_map vm;

  try
    {
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);     
    } 
  catch (const po::error &e)
    {
    std::cout << "Parse error: " << std::endl;
    std::cout << config << std::endl;
    return EXIT_FAILURE;
    }

  // End option reading configuration

  // Display help if asked or program improperly called
  if(vm.count("help") || !vm.count("fiber-file"))
    {
    std::cout << config << std::endl;
    
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      return EXIT_FAILURE;
    }

  const bool VERBOSE = vm.count("verbose");

  // Reader fiber bundle
  GroupType::Pointer group = readFiberFile(vm["fiber-file"].as<std::string>());

  DeformationImageType::Pointer deformationfield(NULL);
  if(vm.count("h-field"))
    deformationfield = readDeformationField(vm["h-field"].as<std::string>(), HField);
//  else if(vm.count("displacement-field"))
//    deformationfield = readDeformationField(vm["displacement-field"].as<std::string>(), Displacement);
  else
    deformationfield = NULL;

  typedef itk::VectorLinearInterpolateImageFunction<DeformationImageType, double> DeformationInterpolateType;
  DeformationInterpolateType::Pointer definterp(NULL);
  if(deformationfield)
    {
    definterp = DeformationInterpolateType::New();
    definterp->SetInputImage(deformationfield);
    }

  // Setup new fiber bundle group
  GroupType::Pointer newgroup = GroupType::New();
  newgroup->SetId(0);

  ChildrenListType* children = group->GetChildren(0);

  if(VERBOSE)
    std::cout << "Getting spacing" << std::endl;

  // Iterate over all tubes in the group
  const double* spacing = dynamic_cast<DTITubeType*>(children->front().GetPointer())->GetSpacing();
  newgroup->SetSpacing(spacing);

  const itk::Vector<double, 3> sooffset = group->GetObjectToParentTransform()->GetOffset();
  newgroup->GetObjectToParentTransform()->SetOffset(sooffset.GetDataPointer());

  // Setup tensor file if available
  typedef itk::DiffusionTensor3D<double> TensorPixelType;
  typedef itk::Image<itk::Vector<double, 6>, 3> TensorImageType;
  typedef itk::Image<itk::Vector<double, 6>, 3> VectorImageType;

  typedef itk::ImageFileReader<TensorImageType> TensorImageReader;
  typedef itk::VectorBSplineInterpolateImageFunction<TensorImageType, double> TensorInterpolateType;
  TensorImageReader::Pointer tensorreader = NULL;
  TensorInterpolateType::Pointer tensorinterp = NULL;
  
  if(vm.count("tensor-volume"))
    {
    tensorreader = TensorImageReader::New();
    tensorinterp = TensorInterpolateType::New();

    tensorreader->SetFileName(vm["tensor-volume"].as<std::string>().c_str());
    try
      {
      tensorreader->Update();
      tensorinterp->SetInputImage(tensorreader->GetOutput());
      }
    catch(itk::ExceptionObject exp)
      {
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
      }
    }

  if(VERBOSE)
    std::cout << "Starting Loop" << std::endl;

  ChildrenListType::iterator it;
  unsigned int id = 1;

  // For each fiber
  for(it = children->begin(); it != children->end(); it++)
    {
    DTIPointListType pointlist = dynamic_cast<DTITubeType*>((*it).GetPointer())->GetPoints();
    DTITubeType::Pointer newtube = DTITubeType::New();
    DTIPointListType newpoints;
    
    DTIPointListType::iterator pit;

    // For each point alogng thje fiber
    for(pit = pointlist.begin(); pit != pointlist.end(); ++pit)
      {
      DTIPointType newpoint;
      typedef DTIPointType::PointType PointType;
      PointType p = pit->GetPosition();
      typedef DeformationInterpolateType::ContinuousIndexType ContinuousIndexType;
      ContinuousIndexType ci;
      for(unsigned int i =0; i < 3; i++)
        ci[i] = p[i];

      if(deformationfield)
        {
        
        DeformationPixelType warp(definterp->EvaluateAtContinuousIndex(ci).GetDataPointer());
        
        if(!vm.count("no-warp"))
          {
          for(unsigned int i =0; i < 3; i++)
            p[i] = p[i] + warp[i] / spacing[i];
          }
        }
      
      newpoint.SetPosition(p);
      newpoint.SetRadius(.4);
      newpoint.SetRed(0.0);
      newpoint.SetGreen(1.0);
      newpoint.SetBlue(0.0);

      // Attribute tensor data if provided
      if(vm.count("tensor-volume"))
        {
        itk::DiffusionTensor3D<double> tensor(tensorinterp->EvaluateAtContinuousIndex(ci).GetDataPointer());

        // TODO: Change SpatialObject interface to accept DiffusionTensor3D
        float sotensor[6];
        for(unsigned int i = 0; i < 6; ++i)
          sotensor[i] = tensor[i];

        newpoint.SetTensorMatrix(sotensor);

        typedef itk::DiffusionTensor3D<double>::EigenValuesArrayType EigenValuesType;
        EigenValuesType eigenvalues;
        tensor.ComputeEigenValues(eigenvalues);
        
        newpoint.AddField(itk::DTITubeSpatialObjectPoint<3>::FA, tensor.GetFractionalAnisotropy());
        newpoint.AddField("md", tensor.GetTrace()/3);
        newpoint.AddField("l1", eigenvalues[0]);
        newpoint.AddField("l2", eigenvalues[1]);
        newpoint.AddField("l3", eigenvalues[2]);
        }
      else
        {
        newpoint.SetTensorMatrix(pit->GetTensorMatrix());
        }

      newpoints.push_back(newpoint);
      }

    newtube->SetSpacing(spacing);
    newtube->SetId(id++);
    newtube->SetPoints(newpoints);
    newgroup->AddSpatialObject(newtube);
    }

  newgroup->ComputeObjectToWorldTransform();

  if(VERBOSE)
    std::cout << "Ending Loop" << std::endl;

  if(VERBOSE)
    std::cout << "Output: " << vm["fiber-output"].as<std::string>() << std::endl;

  writeFiberFile(vm["fiber-output"].as<std::string>(), newgroup);

  delete children;
}
