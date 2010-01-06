/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.3 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <string>
#include <iostream>
#include <cstdlib>

#include <itkInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>

#include <itkAffineTransform.h>

#include <itkResampleImageFilter.h>
#include <itkWarpImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>

#include "deformationfieldio.h"
#include "dtitypes.h"
#include "scalartransformCLP.h"
#if 0
// Validates the interpolation type option string to the the allowed
// values for interpolation methods.  Currently nearestneighbor,
// linear, or cubic.
void validate(boost::any& v,
              const std::vector<std::string>& values,
              InterpolationType* target_type,
              int)
{
  using namespace boost::program_options;;
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
#endif

typedef itk::InterpolateImageFunction<IntImageType, double> InterpolatorType;

InterpolatorType::Pointer createInterpolater(InterpolationType interp)
{
  switch(interp)
    {
    case NearestNeighbor:
      return itk::NearestNeighborInterpolateImageFunction<IntImageType, double>::New().GetPointer();
    case Linear:
      return itk::LinearInterpolateImageFunction<IntImageType, double>::New().GetPointer();
    case Cubic:
      return itk::BSplineInterpolateImageFunction<IntImageType, double>::New().GetPointer();
    default:
      throw itk::ExceptionObject("Invalid interpolation type");
    }
  return NULL;
}

int main(int argc, char* argv[])
{
#if 0
  namespace po = boost::program_options;

  // Read program options/configuration
  po::options_description config("Usage: scalartransform [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v", "Verbose output")

    ("input-image,i", po::value<std::string>(), "Image to transform")
    ("output-image,o", po::value<std::string>(), "The transformed result of the moving image")
    ("transformation,t", po::value<std::string>(), "Output file for transformation parameters")
    ("invert", po::value<bool>()->default_value(false), "Invert transform before applying (default: false)")
    ("deformation,d", po::value<std::string>(), "Deformation Field")
    ("h-field", "The deformation is an h-field")
    ("interpolator", po::value<InterpolationType>()->default_value(Linear, "linear"), "Interpolation type (neareastneighbor, linear, cubic")
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
    std::cout << config << std::endl;
    return EXIT_FAILURE;
    }

  if(vm.count("help") || !vm.count("input-image") || !vm.count("output-image")
     || (!vm.count("transformation") && !vm.count("deformation")))
    {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      {
      std::cerr << "The input, output, and transformation must be specified." << std::endl;
      return EXIT_FAILURE;
      }
    }
#endif
  PARSE_ARGS;
  if(inputImage == "" || outputImage == "" ||
     transformation == "" || deformation == "")
    {
    std::cerr << "The input, output, and transformation must be specified." << std::endl;
    return EXIT_FAILURE;
    }
  typedef itk::ImageFileReader<IntImageType> ImageReader;
  
  ImageReader::Pointer reader = ImageReader::New();
  
  reader->SetFileName( inputImage );

  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  InterpolationType interpType =
    (interpolation == "linear" ? Linear :
     (interpolation == "nearestneightbor" ? NearestNeighbor :
      Cubic));
  IntImageType::Pointer result = NULL;
  InterpolatorType::Pointer interp = createInterpolater(interpType);
  if(transformation != "")
  {
    typedef itk::TransformFileReader TransformReader;
    TransformReader::Pointer treader = TransformReader::New();

    treader->SetFileName(transformation);

    try
    {
      treader->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
    }

    typedef itk::ResampleImageFilter<IntImageType, IntImageType, double> ResampleFilter;
    ResampleFilter::Pointer  resampler = ResampleFilter::New();
    resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
    resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
    resampler->SetInterpolator( interp );

    resampler->SetInput( reader->GetOutput() );

    typedef itk::AffineTransform<double, 3> TransformType;
    TransformType::Pointer transform = dynamic_cast<TransformType*>( treader->GetTransformList()->front().GetPointer() );
    assert(!transform.IsNull());
    resampler->SetTransform( transform );
    resampler->Update();
    result = resampler->GetOutput();
  }
  else if(deformation != "")
  {
    DeformationImageType::Pointer defimage = DeformationImageType::New();
    defimage = readDeformationField(deformation,
                                    hField ? HField : Displacement);

    typedef itk::WarpImageFilter<IntImageType, IntImageType, DeformationImageType> WarpFilter;
    WarpFilter::Pointer warpresampler = WarpFilter::New();
    warpresampler->SetInterpolator(interp);
    warpresampler->SetEdgePaddingValue(0);
    warpresampler->SetInput(reader->GetOutput());
    warpresampler->SetDeformationField(defimage);
    warpresampler->SetOutputSpacing( reader->GetOutput()->GetSpacing());
    warpresampler->SetOutputOrigin( reader->GetOutput()->GetOrigin());
    warpresampler->Update();

    result = warpresampler->GetOutput();
                                       
  }
  else
  {
    std::cerr << "Unknown transformation type" << std::endl;
    return EXIT_FAILURE;
  }
     
  typedef itk::ImageFileWriter<IntImageType > ImageWriter;
  ImageWriter::Pointer writer = ImageWriter::New();
  writer->UseCompressionOn();
  writer->SetFileName(outputImage);
  writer->SetInput(result);

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
