/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $
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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkNthElementImageAdaptor.h>
#include <itkAbsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkVersion.h>
#include "maxcurvatureCLP.h"

enum CurvatureType { MaxEigenvalue, SmoothNormalized, RawNormalized, UnPossible };

int main(int argc, char* argv[])
{
  PARSE_ARGS;

  // Display help if asked or program improperly called
  if( image == "" )
    {
    std::cout << "Version: $Date: 2009-01-09 15:39:51 $ $Revision: 1.5 $" << std::endl;
    std::cout << ITK_SOURCE_VERSION << std::endl;
    return EXIT_FAILURE;
    }

  if( output == "" )
    {
    std::cerr << "maxcurvature: Must specify output file" << std::endl;
    return EXIT_FAILURE;
    }

  typedef unsigned short PixelType;
  typedef double         FloatPixelType;
  const int DIM = 3;
  typedef itk::SymmetricSecondRankTensor<double, DIM> HessianPixelType;
  typedef itk::Vector<double, DIM>                    EigenValuePixelType;

  typedef itk::Image<PixelType, DIM>      ImageType;
  typedef itk::Image<FloatPixelType, DIM> FloatImageType;

  typedef itk::Image<HessianPixelType, DIM>    HessianImageType;
  typedef itk::Image<EigenValuePixelType, DIM> EigenValueImageType;

  typedef itk::ImageFileReader<ImageType> FileReaderType;

  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(image);

  // sigma set by PARSE_ARGS
  //  double sigma = vm["sigma"].as<double>();
  typedef itk::HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;
  HessianFilterType::Pointer hessian = HessianFilterType::New();
  hessian->SetInput(reader->GetOutput() );
  hessian->SetSigma(sigma);

  typedef itk::SymmetricEigenAnalysisImageFilter<HessianImageType, EigenValueImageType> EigenAnalysisFilterType;
  EigenAnalysisFilterType::Pointer eigen = EigenAnalysisFilterType::New();
  eigen->SetInput(hessian->GetOutput() );
  eigen->OrderEigenValuesBy(EigenAnalysisFilterType::FunctorType::OrderByValue);
  eigen->SetDimension(3);

  eigen->Update();
  typedef itk::NthElementImageAdaptor<EigenValueImageType, FloatPixelType> ElementSelectAdaptorType;
  ElementSelectAdaptorType::Pointer elementSelect = ElementSelectAdaptorType::New();
  elementSelect->SetImage(eigen->GetOutput() );
  elementSelect->SelectNthElement(0);

  typedef itk::CastImageFilter<ElementSelectAdaptorType, FloatImageType> CastFilter2Type;
  CastFilter2Type::Pointer cast2 = CastFilter2Type::New();
  cast2->SetInput(elementSelect);
  cast2->Update();

  typedef itk::ShiftScaleImageFilter<FloatImageType, FloatImageType> ScaleImageType;
  ScaleImageType::Pointer scale = ScaleImageType::New();
  scale->SetShift(0.0);

  scale->SetScale(-10.0);
  scale->SetInput(cast2->GetOutput() );

  typedef itk::IntensityWindowingImageFilter<FloatImageType> WindowFilterType;
  WindowFilterType::Pointer window = WindowFilterType::New();
  window->SetInput(scale->GetOutput() );
  window->SetWindowMinimum(0);
  window->SetWindowMaximum(32768);
  window->SetOutputMinimum(0);
  window->SetOutputMaximum(32768);

  typedef itk::CastImageFilter<FloatImageType, ImageType> CastFilterType;
  CastFilterType::Pointer cast = CastFilterType::New();
  cast->SetInput(window->GetOutput() );

  typedef itk::ImageFileWriter<ImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetUseCompression(true);
  writer->SetInput(cast->GetOutput() );
  writer->SetFileName(output);

  try
    {
    eigen->Update();
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    return -1;
    }

  return 0;
}
