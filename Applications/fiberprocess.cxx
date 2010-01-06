/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/08/03 17:36:42 $
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
#include <fstream>

// ITK includes
#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTensorLinearInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkVersion.h>

//#include "FiberCalculator.h"
#include "deformationfieldoperations.h"
#include "fiberio.h"
#include "dtitypes.h"
#include "fiberprocessCLP.h"

int main(int argc, char* argv[])
{
#if 0
  // Read program options/configuration
  po::options_description config("Usage: fiberprocess input-fiber [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v", "produces verbose output")

    ("fiber-output,o", po::value<std::string>(), "Output fiber file.  May be warped or updated with new data depending on other options used.")
    ("h-field,H", po::value<std::string>(), "HField for warp and statistics lookup.  If this option is used tensor-volume must also be specified.")
    ("displacement-field", po::value<std::string>(), "Displacement field for warp and statistics lookup.  If this option is used tensor-volume must also be specified.")
    ("no-warp,n", "Do not warp the geometry of the tensors only obtain the new statistics")
    ("index-space", "Use index-space for fiber output coordinates, otherwise us world space for fiber output coordinates (from tensor file)")
    ("tensor-volume,T", po::value<std::string>(), "Interpolate tensor values from the given field")

    // ******  TODO **********
    ("voxelize,V", po::value<std::string>(),"Voxelize fiber into a label map (the labelmap filename is the argument of -V).  The tensor file must be specified using -T for information about the size, origin, spacing of the image.")
    ("voxelize-count-fibers", "Count number of fibers per-voxel instead of just setting to 1")
    ("voxel-label,l", po::value<ScalarPixelType>()->default_value(1),"Label for voxelized fiber")
    ("fiberRadius,R",po::value<double>(),"set radius of all fibers to this value")

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
    {
      std::cout << "Version $Revision: 1.8 $ "<< std::endl;
      std::cout << ITK_SOURCE_VERSION << std::endl;
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
#endif
  PARSE_ARGS;
  const bool VERBOSE = verbose;

  // Reader fiber bundle
  GroupType::Pointer group = readFiberFile(fiberFile);

  DeformationImageType::Pointer deformationfield(NULL);
  if(hField != "")
    deformationfield = readDeformationField(hField, HField);
  else if(displacementField != "")
    deformationfield = readDeformationField(displacementField, Displacement);
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

  // Get Spacing and offset from group
  double spacing[3];
  for(unsigned int i =0; i < 3; i++)
    spacing[i] = (group->GetSpacing())[i];
  newgroup->SetSpacing(spacing);

  itk::Vector<double, 3> sooffset;
  for(unsigned int i = 0; i < 3; i++)
    sooffset[i] = (group->GetObjectToParentTransform()->GetOffset())[i];
  newgroup->GetObjectToParentTransform()->SetOffset(sooffset.GetDataPointer());

  if(VERBOSE)
  {
    std::cout << "Group Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;
    std::cout << "Group Offset: " << sooffset[0] << ", " << sooffset[1]  << ", " << sooffset[2] << std::endl;
  }

  // Setup tensor file if available
  typedef itk::ImageFileReader<TensorImageType> TensorImageReader;
  typedef itk::TensorLinearInterpolateImageFunction<TensorImageType, double> TensorInterpolateType;
  TensorImageReader::Pointer tensorreader = NULL;
  TensorInterpolateType::Pointer tensorinterp = NULL;
  
  if(tensorVolume != "")
  {
    tensorreader = TensorImageReader::New();
    tensorinterp = TensorInterpolateType::New();

    tensorreader->SetFileName(tensorVolume);
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

  // Need to allocate an image to write into for creating
  // the fiber label map
  IntImageType::Pointer labelimage;
  if(voxelize != "")
  {
    if(tensorVolume == "")
    {
      std::cerr << "Must specify tensor file to copy image metadata for fiber voxelize." << std::endl;
      return EXIT_FAILURE;
    }
    tensorreader->GetOutput();
    labelimage = IntImageType::New();
    labelimage->SetSpacing(tensorreader->GetOutput()->GetSpacing());
    labelimage->SetOrigin(tensorreader->GetOutput()->GetOrigin());
    labelimage->SetDirection(tensorreader->GetOutput()->GetDirection());
    labelimage->SetRegions(tensorreader->GetOutput()->GetLargestPossibleRegion());
    labelimage->Allocate();
    labelimage->FillBuffer(0);
  }

  double newspacing[3];
  for(unsigned int i =0; i < 3; i++)
    newspacing[i] = spacing[i];
  
  if (tensorVolume != "" || voxelize != "")
  {
    for(unsigned int i =0; i < 3; i++)
      newspacing[i] = 1;
    // fiber coordinates are newly in mm
  }

  // For each fiber
  for(it = children->begin(); it != children->end(); it++)
  {
    // Get Spacing and offset from fiber
    const double* fiberspacing = (*it).GetPointer()->GetSpacing();
    const itk::Vector<double, 3> fibersooffset = dynamic_cast<DTITubeType*>((*it).GetPointer())->GetObjectToParentTransform()->GetOffset();
    
    // TODO: mode for selecting which to use
    if ((fiberspacing[0] != spacing[0]) || (fiberspacing[1] != spacing[1]) || (fiberspacing[1] != spacing[1]) || 
	(sooffset[0] != fibersooffset[0]) || (sooffset[1] != fibersooffset[1]) || (sooffset[1] != fibersooffset[1]))
    {
      std::cout << "Fiber and group spacing/offset is not equal using the one from the fiber group." << std::endl;
      if(VERBOSE)
      {
	std::cout << "Fiber Spacing: " << fiberspacing[0] << ", " << fiberspacing[1] << ", " << fiberspacing[2] << std::endl;
	std::cout << "Fiber Offset: " << fibersooffset[0] << ", " << fibersooffset[1]  << ", " << fibersooffset[2] << std::endl;
      }
      for(unsigned int i =0; i < 3; i++) 
      {
	spacing[i] = fiberspacing[i];
	sooffset[i] = fibersooffset[i];
      }
      newgroup->SetSpacing(spacing);
      newgroup->GetObjectToParentTransform()->SetOffset(sooffset.GetDataPointer());
    }

    DTIPointListType pointlist = dynamic_cast<DTITubeType*>((*it).GetPointer())->GetPoints();
    DTITubeType::Pointer newtube = DTITubeType::New();
    DTIPointListType newpoints;
    
    DTIPointListType::iterator pit;

    // For each point along the fiber
    for(pit = pointlist.begin(); pit != pointlist.end(); ++pit)
    {
      DTIPointType newpoint;
      typedef DTIPointType::PointType PointType;

      // p is not really a point its a continuous index
      const PointType p = pit->GetPosition();
      itk::Point<double, 3> pt, newpt, trypt; 
      // pt is Point in world coordinates (without any reorientation, i.e. disregarding transform matrix)
      pt[0] = p[0] * spacing[0] + sooffset[0];
      pt[1] = p[1] * spacing[1] + sooffset[1];
      pt[2] = p[2] * spacing[2] + sooffset[2];

      typedef DeformationInterpolateType::ContinuousIndexType ContinuousIndexType;
      ContinuousIndexType ci, origci;
      for(unsigned int i =0; i < 3; i++)
        origci[i] = ci[i] = p[i];

      newpoint.SetPosition(ci);

      if (tensorVolume != "")
      {
	// get index in image
	tensorreader->GetOutput()->TransformPhysicalPointToContinuousIndex(pt, ci);
	for(unsigned int i =0; i < 3; i++)
	  origci[i] = ci[i];
	if(deformationfield)
	{
	  DeformationPixelType warp(definterp->EvaluateAtContinuousIndex(ci).GetDataPointer());
        
	  // Get Spacing from tensorfile
	  TensorImageType::SpacingType tensorspacing = tensorreader->GetOutput()->GetSpacing();

	  for(unsigned int i =0; i < 3; i++)
	    ci[i] = ci[i] + warp[i] / tensorspacing[i];

	  if(indexSpace)
	  {
          if(noWarp)
            {
            newpoint.SetPosition(origci);
            }
          else
            {
            newpoint.SetPosition(ci);
            }
	  }
	  else
	  {
	    tensorreader->GetOutput()->TransformContinuousIndexToPhysicalPoint(ci, newpt);
	    // this uses full physical space info from tensor image file
	    if(noWarp)
              {
	      newpoint.SetPosition(pt);
              }
	    else
              {
	      newpoint.SetPosition(newpt);
              }
	  }
	  //	  double  newFiberSpacing[3];
	  //	  for(unsigned int i =0; i < 3; i++)
	  //	    newFiberSpacing[i] =  tensorspacing[i];	  
	  //	  newgroup->SetSpacing(newFiberSpacing);
	} 
      }

      if(voxelize != "")
      {
        itk::Index<3> ind;
        labelimage->TransformPhysicalPointToContinuousIndex(pt, ci);
        ind[0] = static_cast<long int>(round(ci[0]));
        ind[1] = static_cast<long int>(round(ci[1]));
        ind[2] = static_cast<long int>(round(ci[2]));
        if(!labelimage->GetLargestPossibleRegion().IsInside(ind))
        {
          std::cerr << "Error index: " << ind << " not in image"  << std::endl;
          return EXIT_FAILURE;
        }
        if(voxelizeCountFibers)
          {
          labelimage->SetPixel(ind, labelimage->GetPixel(ind) + 1);
          }
        else
          {
          labelimage->SetPixel(ind, voxelLabel);
          }
      }
#if 0
      if (vm.count("fiberRadius")) {
	newpoint.SetRadius(vm["fiberRadius"].as<double>());
      } 
      else
      {
	newpoint.SetRadius(0.4);
      }
#else
      newpoint.SetRadius(fiberRadius);
#endif
      newpoint.SetRed(0.0);
      newpoint.SetGreen(1.0);
      newpoint.SetBlue(0.0);

      // Attribute tensor data if provided
      if(tensorVolume != "" && fiberOutput != "")
      { 
	//std::cout << "ci" << ci[0] << "," << ci[1] << "," << ci[2] << std::endl;
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
        newpoint.AddField("fro", sqrt(tensor[0]*tensor[0] +
                                      2*tensor[1]*tensor[1] +
                                      2*tensor[2]*tensor[2] +
                                      tensor[3]*tensor[3] +
                                      2*tensor[4]*tensor[4] +
                                      tensor[5]*tensor[5]));
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

    newtube->SetSpacing(newspacing);
    newtube->SetId(id++);
    newtube->SetPoints(newpoints);
    newgroup->AddSpatialObject(newtube);
  }

  newgroup->ComputeObjectToWorldTransform();

  if(VERBOSE)
    std::cout << "Ending Loop" << std::endl;

  if(VERBOSE && fiberOutput != "")
    std::cout << "Output: " << fiberOutput << std::endl;

  if(fiberOutput != "")
    writeFiberFile(fiberOutput, newgroup);

  if(voxelize != "")
  {
    typedef itk::ImageFileWriter<IntImageType> LabelWriter;
    LabelWriter::Pointer writer = LabelWriter::New();
    writer->SetInput(labelimage);
    writer->SetFileName(voxelize);
    writer->UseCompressionOn();
    try
    {
      writer->Update();
    }
    catch(itk::ExceptionObject & e)
    {
      std::cerr << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  }

  delete children;
  return EXIT_SUCCESS;
}
