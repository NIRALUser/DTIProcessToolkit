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

typedef itk::InterpolateImageFunction<IntImageType, double> InterpolatorType;

InterpolatorType::Pointer createInterpolater(InterpolationType interp)
{
  switch( interp )
    {
    case NearestNeighbor:
      {
      return itk::NearestNeighborInterpolateImageFunction<IntImageType, double>::New().GetPointer();
      }
    case Linear:
      {
      return itk::LinearInterpolateImageFunction<IntImageType, double>::New().GetPointer();
      }
    case Cubic:
      {
      return itk::BSplineInterpolateImageFunction<IntImageType, double>::New().GetPointer();
      }
    default:
      throw itk::ExceptionObject("Invalid interpolation type");
    }
  return ITK_NULLPTR;
}

int main(int argc, char* argv[])
{
  PARSE_ARGS;
  if( inputImage == "" || outputImage == "" ||
      transformation == "" || deformation == "" )
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
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  InterpolationType interpType =
    (interpolation == "linear" ? Linear :
     (interpolation == "nearestneightbor" ? NearestNeighbor :
      Cubic) );
  IntImageType::Pointer     result = ITK_NULLPTR;
  InterpolatorType::Pointer interp = createInterpolater(interpType);
  if( transformation != "" )
    {
    typedef itk::TransformFileReader TransformReader;
    TransformReader::Pointer treader = TransformReader::New();

    treader->SetFileName(transformation);

    try
      {
      treader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }

    typedef itk::ResampleImageFilter<IntImageType, IntImageType, double> ResampleFilter;
    ResampleFilter::Pointer resampler = ResampleFilter::New();
    resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
    resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
    resampler->SetInterpolator( interp );

    resampler->SetInput( reader->GetOutput() );

    typedef itk::AffineTransform<double, 3> TransformType;
    TransformType::Pointer transform =
      dynamic_cast<TransformType *>( treader->GetTransformList()->front().GetPointer() );
    assert(!transform.IsNull() );
    resampler->SetTransform( transform );
    resampler->Update();
    result = resampler->GetOutput();
    }
  else if( deformation != "" )
    {
    DeformationImageType::Pointer defimage = DeformationImageType::New();
    defimage = readDeformationField(deformation,
                                    hField ? HField : Displacement);

    typedef itk::WarpImageFilter<IntImageType, IntImageType, DeformationImageType> WarpFilter;
    WarpFilter::Pointer warpresampler = WarpFilter::New();
    warpresampler->SetInterpolator(interp);
    warpresampler->SetEdgePaddingValue(0);
    warpresampler->SetInput(reader->GetOutput() );
#if ITK_VERSION_MAJOR < 4
    warpresampler->SetDeformationField(defimage);
#else
    warpresampler->SetDisplacementField(defimage);
#endif
    warpresampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
    warpresampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
    warpresampler->Update();

    result = warpresampler->GetOutput();

    }
  else
    {
    std::cerr << "Unknown transformation type" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<IntImageType> ImageWriter;
  ImageWriter::Pointer writer = ImageWriter::New();
  writer->UseCompressionOn();
  writer->SetFileName(outputImage);
  writer->SetInput(result);

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
