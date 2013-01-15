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

#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>

int main(int argc, char* argv[])
{
  if(argc != 5)
  {
    std::cerr << "Usage: " << argv[0] << " type radius input output" << std::endl;
    return EXIT_FAILURE;
  }

  const std::string type = argv[1];
  const unsigned int radius = atoi(argv[2]);
  const std::string input = argv[3];
  const std::string output = argv[4];
 
  typedef unsigned short Pixel;
  typedef itk::Image<Pixel, 3> Image;
  typedef itk::BinaryBallStructuringElement<Pixel, 3> StructuringElement;

  itk::ImageToImageFilter<Image, Image>::Pointer filter;
  StructuringElement ball;
  ball.SetRadius(radius);
  ball.CreateStructuringElement();

  if(type == "open")
  {
    typedef itk::BinaryMorphologicalOpeningImageFilter<Image, Image, StructuringElement> OpeningFilter;
    OpeningFilter::Pointer openingfilter = OpeningFilter::New();
    openingfilter->SetKernel(ball);
    openingfilter->SetForegroundValue(1);
    filter = openingfilter;
  }
  else if(type == "close")
  {
    typedef itk::BinaryMorphologicalClosingImageFilter<Image, Image, StructuringElement> ClosingFilter;
    ClosingFilter::Pointer closingfilter = ClosingFilter::New();
    closingfilter->SetKernel(ball);
    closingfilter->SetForegroundValue(1);
    filter = closingfilter;
  }
  else if(type == "dilate")
  {
    typedef itk::BinaryDilateImageFilter<Image, Image, StructuringElement> DilateFilter;
    DilateFilter::Pointer dilatefilter = DilateFilter::New();
    dilatefilter->SetKernel(ball);
    dilatefilter->SetForegroundValue(1);
    filter = dilatefilter;
  }
  else if(type == "erode")
  {
    typedef itk::BinaryErodeImageFilter<Image, Image, StructuringElement> ErodeFilter;
    ErodeFilter::Pointer erodefilter = ErodeFilter::New();
    erodefilter->SetKernel(ball);
    erodefilter->SetForegroundValue(1);
    filter = erodefilter;
  }
  else
  {
    std::cerr << "Invalid morpholigical operation" << std::endl;
  }
    
  typedef itk::ImageFileReader<Image> ImageReader;
  typedef itk::ImageFileWriter<Image> ImageWriter;

  ImageReader::Pointer reader = ImageReader::New();
  reader->SetFileName(input);
  ImageWriter::Pointer writer = ImageWriter::New();
  writer->SetFileName(output);
  writer->UseCompressionOn();

  filter->SetInput(reader->GetOutput());
  writer->SetInput(filter->GetOutput());

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject e)
  {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
