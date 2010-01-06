/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.7 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>
#include <exception>
#include <cassert>

#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkImageToVectorImageFilter.h>
#include <itkNthElementImageAdaptor.h>

#include "itkVectorBSplineInterpolateImageFunction.h"
#include "transforms.h"

int main(int argc, char* argv[])
{
  // This software reads a vectorized .nrrd file performs a
  // transformation and updates the embedded gradient strings

  typedef double RealType;
  typedef double TransformRealType;
  typedef itk::AffineTransform<TransformRealType,3> AffineTransformType;

  const unsigned int DIM = 3;
  typedef unsigned short DWIPixelType;
  typedef itk::VectorImage<DWIPixelType, DIM> VectorImageType;
  typedef itk::Image<DWIPixelType, DIM> ComponentImageType;

  if(argc != 4)
    {
    std::cerr << "Usage: " << argv[0] << " infile outfile transform" << std::endl;
    std::cerr << "This program transform a vector image and updates the gradient"
              << " directions if they are specified using the NAMIC convention for DWI data";

    return EXIT_FAILURE;
    }

  const std::string infile = argv[1];
  const std::string outfile = argv[2];
  const std::string transformfile = argv[3];

  typedef itk::ImageFileReader<VectorImageType> ImageFileReaderType;
  ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
  reader->SetFileName(infile);
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e <<std::endl;
    return EXIT_FAILURE;
    }

  VectorImageType::Pointer dwimg = reader->GetOutput();

  AffineTransformType::Pointer transform = NULL;
  vnl_matrix<TransformRealType> R(3,3,0);
  if (transformfile.rfind(".dof") != std::string::npos)
    {
    // Use RView transform file reader
    RViewTransform<TransformRealType> dof(readDOFFile<TransformRealType>(transformfile));
    // image transform
    transform = createITKAffine(dof,
                                dwimg->GetLargestPossibleRegion().GetSize(),
                                dwimg->GetSpacing(),
                                dwimg->GetOrigin());

    // gradient transform
    // g' = R g
    R = vnl_matrix_inverse<TransformRealType>(transform->GetMatrix().GetVnlMatrix());

    }
  else
    {
    // Assume ITK transform
    std::cerr << "ITK transform not implemented" << std::endl;
    return EXIT_FAILURE;
    }

  // Transform components of dwi image
  // TODO: this really isnt that awesome or memory efficient a way of
  // doing things.
  const unsigned int vectorlength = reader->GetOutput()->GetVectorLength();
  assert(vectorlength > 6);
  
  typedef itk::NthElementImageAdaptor<VectorImageType, unsigned int> VectorElementAdaptorType;
  VectorElementAdaptorType::Pointer adaptor = VectorElementAdaptorType::New();
  adaptor->SetImage(reader->GetOutput());
  
  typedef itk::ImageToVectorImageFilter<ComponentImageType> VectorComposeFilterType;
  VectorComposeFilterType::Pointer vcompose = VectorComposeFilterType::New();

  typedef itk::ResampleImageFilter<VectorElementAdaptorType,ComponentImageType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();

  for(unsigned int i = 0; i < vectorlength; i++)
    {
    // Transform the ith component of the vector image
    adaptor->SelectNthElement(i);
    adaptor->Update();

    resampler->SetInput(adaptor);
    resampler->SetTransform(transform);
    resampler->SetOutputSpacing(reader->GetOutput()->GetSpacing());
    resampler->SetOutputOrigin(reader->GetOutput()->GetOrigin());
    resampler->SetSize(reader->GetOutput()->GetLargestPossibleRegion().GetSize());
    resampler->Update();

    vcompose->SetInput(i,resampler->GetOutput());
    }
  vcompose->Update();

  // Transform gradient directions

  // write output
  typedef itk::ImageFileWriter<VectorImageType> ImageFileWriterType;
  ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
  writer->SetInput(vcompose->GetOutput());
  writer->SetFileName(outfile);
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
