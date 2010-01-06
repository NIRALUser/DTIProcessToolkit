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
#include <iostream>
#include <fstream>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkChangeLabelImageFilter.h>

int main(int argc, char* argv[])
{
  const char* infile = argv[1];
  const char* outfile = argv[2];

  typedef unsigned short PixelType;
  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageFileWriterType;

  typedef itk::ChangeLabelImageFilter<ImageType, ImageType> ChangeLabelType;

  ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
  reader->SetFileName(infile);

  ChangeLabelType::Pointer labelchanger = ChangeLabelType::New();
  labelchanger->SetInput(reader->GetOutput());

  for(int i = 3; i  < argc;  i += 2)
    labelchanger->SetChange(atoi(argv[i]),atoi(argv[i+1]));

  ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
  writer->SetUseCompression(true);
  writer->SetInput(labelchanger->GetOutput());
  writer->SetFileName(outfile);
  writer->Update();

  return EXIT_SUCCESS;
}
