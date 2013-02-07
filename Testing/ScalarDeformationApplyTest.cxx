#include <iostream>
#include <fstream>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkResampleImageFilter.h>
#include <itkWarpImageFilter.h>

#include "itkDeformationFieldFromTransform.h"

#include "transforms.h"

int ScalarDeformationApplyTest(int argc, char* argv[])
{
  // 1: input file
  // 2: transform
  // 3: transform image
  // 4: warped image
  if(argc < 5)
    {
    std::cerr << argv[0] << ": input transform transformImage warpedImage" << std::endl;
    return 1;
    }
  typedef double TransformRealType;
  typedef unsigned short PixelType;
  typedef itk::Image<PixelType,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  RViewTransform<TransformRealType> dof = readDOFFile<TransformRealType>(argv[2]);
  typedef itk::AffineTransform<TransformRealType,3> AffineTransformType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  ImageType::Pointer image = reader->GetOutput();

  typedef itk::Vector<double,3> DeformationPixelType;
  typedef itk::Image<DeformationPixelType,3> DeformationImageType;

  typedef DeformationImageType::SizeType ImageSizeType;
  typedef DeformationImageType::SpacingType ImageSpacingType;
  typedef itk::DeformationFieldFromTransform<DeformationImageType,double> DeformationSourceType;
  
  DeformationSourceType::Pointer defgen = DeformationSourceType::New();

  defgen->SetOutputRegion(image->GetLargestPossibleRegion());

  AffineTransformType::Pointer transform = createITKAffine(dof,
                                                           image->GetLargestPossibleRegion().GetSize(),
                                                           image->GetSpacing(),
                                                           image->GetOrigin());
  defgen->SetTransform(transform);
  defgen->SetOutputSpacing(image->GetSpacing());
  defgen->SetOutputOrigin(image->GetOrigin());
  defgen->Update();

  // Apply transform using original transform
  typedef itk::ResampleImageFilter<ImageType,ImageType,double> ResampleImageType;
  ResampleImageType::Pointer resample = ResampleImageType::New();
  resample->SetInput(image);
  resample->SetTransform(transform);
  resample->SetOutputParametersFromImage(image);
  resample->Update();

  typedef itk::ImageFileWriter<ImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetInput(resample->GetOutput());
  writer->SetFileName(argv[3]);
  writer->Update();

  // Apply transform using 
  typedef itk::WarpImageFilter<ImageType,ImageType,DeformationImageType> WarpImageType;
  WarpImageType::Pointer warp = WarpImageType::New();
  warp->SetInput(image);
  warp->SetOutputSpacing(image->GetSpacing());
  warp->SetOutputOrigin(image->GetOrigin());
  warp->SetDeformationField(defgen->GetOutput());
  warp->Update();

  writer->SetInput(warp->GetOutput());
  writer->SetFileName(argv[4]);
  writer->Update();

  return 0;
}
