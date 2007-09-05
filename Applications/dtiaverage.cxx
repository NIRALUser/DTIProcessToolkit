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
#include <iostream>
#include <string>
#include <vector>

#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVector.h>
#include <itkExtractImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkImageDuplicator.h>

#include "itkLogEuclideanTensorImageFilter.h"
#include "itkExpEuclideanTensorImageFilter.h"

class PixelDivider
{
public:
  // Divide by zero so we get warned if we didnt replace the default functor
  PixelDivider() : m_Denominator(0) {}
  PixelDivider(double denominator) : m_Denominator(denominator) {}

  itk::Vector<double, 6> operator()(const itk::Vector<double,6> &numerator) { return numerator / m_Denominator; }

  bool operator !=(const PixelDivider & rhs) { return true; }

  double m_Denominator;
};

int main(int argc, char* argv[])
{
  if(argc < 3)
    {
    std::cout << "Usage: " << argv[0] << " outputfile inputfiles..." << std::endl;
    }

  std::string ofile = argv[1];
  std::vector<std::string> sources;
  for(int i = 2; i < argc; ++i)
    {
    sources.push_back(argv[i]);
    }

  typedef double RealType;
  typedef itk::DiffusionTensor3D<RealType> TensorPixelType;
  typedef itk::Image<TensorPixelType, 3> TensorImageType;
  typedef itk::ImageFileReader<TensorImageType> TensorFileReader;
  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType LogTensorImageType;
  
  typedef itk::AddImageFilter<LogTensorImageType, LogTensorImageType, LogTensorImageType> AddImageFilter;
  typedef itk::UnaryFunctorImageFilter<LogTensorImageType, LogTensorImageType, PixelDivider> DivideImageFilter;

  typedef itk::ImageDuplicator<LogTensorImageType> DuplicateImageFilter;

  int n = sources.size();
  
  TensorFileReader::Pointer reader = TensorFileReader::New();
  AddImageFilter::Pointer adder = AddImageFilter::New();
  LogEuclideanFilter::Pointer logfilt = LogEuclideanFilter::New();
  logfilt->SetNumberOfThreads(1);

  DuplicateImageFilter::Pointer dup = DuplicateImageFilter::New();

  reader->SetFileName(sources[0].c_str());
  reader->Update();
  logfilt->SetInput(reader->GetOutput());
  logfilt->Update();
  dup->SetInputImage(logfilt->GetOutput());
  dup->Update();
  LogTensorImageType::Pointer average = dup->GetOutput();

  for(int i = 1; i < n; ++i)
    {
    reader->SetFileName(sources[i].c_str());
    logfilt->SetInput(reader->GetOutput());
    logfilt->Update();

    adder->SetInput1(average);
    adder->SetInput2(logfilt->GetOutput());
    adder->Update();

    dup->SetInputImage(adder->GetOutput());
    dup->Update();
    average = dup->GetOutput();
    
    }

  DivideImageFilter::Pointer divide = DivideImageFilter::New();
  divide->SetFunctor(PixelDivider(n));
  divide->SetInput(average);
  divide->Update();

//   std::cout << "Using Extract" << std::endl;

//   typedef itk::ExtractImageFilter<LogTileImageType, LogTensorImageType> ExtractFilterType;
//   ExtractFilterType::Pointer extract = ExtractFilterType::New();
//   LogTileImageType::RegionType inputregion = accum->GetOutput()->GetLargestPossibleRegion();
//   LogTileImageType::SizeType size = inputregion.GetSize();
//   LogTileImageType::IndexType ind = inputregion.GetIndex();
//   size[3] = 0;
//   ind[3] = 0;
//   inputregion.SetSize(size);
//   inputregion.SetIndex(ind);

//   extract->SetInput(accum->GetOutput());
//   extract->SetExtractionRegion(inputregion);
  
  typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
  ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
  expf->SetInput(divide->GetOutput());
  expf->SetNumberOfThreads(1);
//  expf->Update();

  typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;
  TensorFileWriterType::Pointer twrit = TensorFileWriterType::New();
  twrit->SetUseCompression(true);
  twrit->SetInput(expf->GetOutput());
  twrit->SetFileName(ofile.c_str());
  twrit->Update();

}

