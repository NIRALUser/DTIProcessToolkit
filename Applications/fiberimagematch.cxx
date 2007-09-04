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
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkAffineTransform.h>

#include "itkHFieldToDeformationFieldImageFilter.h"


namespace po = boost::program_options;


int main(int argc, char* argv[])
{
  // Read program options/configuration
  po::options_description config("Usage: fiberprocess input-fiber output-fiber [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v","produces verbose output")
    ("h-field,H", po::value<std::string>(), "HField for warp")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("fiber-file", po::value<std::string>(), "DTI fiber file")
    ("fiber-output", po::value<std::string>(), "Warped fiber file based on a deformation field.  Must input h as \"h-field\" of transform")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("fiber-file",1);
  p.add("fiber-output",1);

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

  // End option reading configuration

  // Display help if asked or program improperly called
  if(vm.count("help") || !vm.count("fiber-file") ||
    !vm.count("fiber-output") || !vm.count("h-field"))
    {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      return EXIT_FAILURE;
    }

  typedef itk::SpatialObjectReader<3, unsigned char> SpatialObjectReaderType;

  typedef itk::DTITubeSpatialObject<3> DTITubeType;
  typedef DTITubeType::TubePointType DTIPointType;
  typedef DTITubeType::PointListType DTIPointListType;

  typedef itk::GroupSpatialObject<3> GroupType;
  typedef GroupType::ChildrenListType ChildrenListType;

  // Reading spatial object
  SpatialObjectReaderType::Pointer soreader = SpatialObjectReaderType::New();

  soreader->SetFileName(vm["fiber-file"].as<std::string>().c_str());
  soreader->Update();

  std::string warpfile(vm["h-field"].as<std::string>());

  // Reading deformation field
  typedef itk::Vector<double, 3> DeformationPixelType;
  typedef itk::Image<DeformationPixelType, 3> DeformationImageType;
  typedef itk::ImageFileReader<DeformationImageType> DeformationImageReader;
  DeformationImageReader::Pointer defreader = DeformationImageReader::New();
  defreader->SetFileName(warpfile.c_str());

  typedef itk::HFieldToDeformationFieldImageFilter<DeformationImageType> DeformationConvertType;
  DeformationConvertType::Pointer defconv = DeformationConvertType::New();
  defconv->SetInput(defreader->GetOutput());
//  defconv->SetSpacing(timg->GetSpacing());
  defconv->Update();

  typedef itk::LinearInterpolateImageFunction<DeformationImageType,double> DeformationInterpolateType;
  DeformationInterpolateType::Pointer definterp = DeformationInterpolateType::New();
  definterp->SetInputImage(defconv->GetOutput());

  // Creating new spatial object group
  GroupType::Pointer group = soreader->GetGroup();

  GroupType::Pointer newgroup = GroupType::New();
  newgroup->SetId(0);

  ChildrenListType* children = group->GetChildren(0);

  std::cout << "Getting spacing" << std::endl;

  // Iterate over all tubes in the group
  const double* sospacing = dynamic_cast<DTITubeType*>(children->front().GetPointer())->GetSpacing();
  double spacing[3];
  spacing[0] = sospacing[0];
  spacing[1] = sospacing[1];
  spacing[2] = sospacing[2];

  newgroup->SetSpacing(spacing);

  std::cout << "Starting Loop" << std::endl;

  ChildrenListType::iterator it;
  unsigned int id = 1;
  for(it = children->begin(); it != children->end(); it++)
    {
    DTIPointListType pointlist = dynamic_cast<DTITubeType*>((*it).GetPointer())->GetPoints();
    DTITubeType::Pointer newtube = DTITubeType::New();
    DTIPointListType newpoints(pointlist.size());
    
    std::copy(pointlist.begin(),pointlist.end(),newpoints.begin());

    DTIPointListType::iterator pit;

    // For each point in the tube
    for(pit = newpoints.begin(); pit != newpoints.end(); ++pit)
      {
      typedef DTIPointType::PointType PointType;
      PointType p = pit->GetPosition();
      typedef DeformationInterpolateType::ContinuousIndexType ContinuousIndexType;
      ContinuousIndexType ci;
      for(unsigned int i =0; i < 3; i++)
        ci[i] = p[i];

      DeformationPixelType warp = definterp->EvaluateAtContinuousIndex(ci);
      
      for(unsigned int i =0; i < 3; i++)
        p[i] = p[i] + warp[i] / spacing[i];

      pit->SetPosition(p);
      pit->SetRadius(.5);
      pit->SetRed(0.0);
      pit->SetGreen(1.0);
      pit->SetBlue(0.0);
      }
    newtube->SetSpacing(spacing);
    newtube->SetId(id++);
    newtube->SetPoints(newpoints);
    newgroup->AddSpatialObject(newtube);
    }

  std::cout << "Ending Loop" << std::endl;

  typedef itk::SpatialObjectWriter<3> WriterType;
  WriterType::Pointer writer  = WriterType::New();
  writer->SetInput(newgroup);
  writer->SetFileName(vm["fiber-output"].as<std::string>().c_str());
  writer->Update();

  delete children;
}
