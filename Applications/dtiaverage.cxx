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
#include <itkVersion.h>
#include <itkCastImageFilter.h>

#include "itkLogEuclideanTensorImageFilter.h"
#include "itkExpEuclideanTensorImageFilter.h"
#include "dtiaverageCLP.h"

template <class TElementType>
class PixelDivider
{
public:
  // Divide by zero so we get warned if we didnt replace the default functor
  PixelDivider() : m_Denominator(0)
  {
  }

  PixelDivider(double denominator) : m_Denominator(denominator)
  {
  }

  TElementType operator()(const TElementType & numerator)
  {
    return numerator / m_Denominator;
  }

  bool operator !=(const PixelDivider & rhs)
  {
    return this != &rhs;
  }

  double m_Denominator;
};

enum StatisticsType { Euclidean, LogEuclidean, PGA };
int main(int argc, char* argv[])
{
  PARSE_ARGS;

  typedef double                                       RealType;
  typedef itk::DiffusionTensor3D<RealType>             TensorPixelType;
  typedef itk::Image<TensorPixelType, 3>               TensorImageType;
  typedef itk::ImageFileReader<TensorImageType>        TensorFileReader;
  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType          LogTensorImageType;
  typedef LogTensorImageType::PixelType                LogTensorPixelType;

  typedef itk::AddImageFilter<LogTensorImageType, LogTensorImageType,
                              LogTensorImageType>                         AddImageFilter;
  typedef itk::UnaryFunctorImageFilter<LogTensorImageType, LogTensorImageType,
                                       PixelDivider<LogTensorPixelType> > DivideImageFilter;

  typedef itk::ImageDuplicator<LogTensorImageType> DuplicateImageFilter;

  //  const std::vector<std::string> sources = vm["inputs"].as<std::vector<std::string> >();

  const int numberofinputs = inputs.size();
  if( numberofinputs > 0 )
    {
    TensorFileReader::Pointer   reader = TensorFileReader::New();
    AddImageFilter::Pointer     adder = AddImageFilter::New();
    LogEuclideanFilter::Pointer logfilt = LogEuclideanFilter::New();

    DuplicateImageFilter::Pointer dup = DuplicateImageFilter::New();
    std::cout << "lbl1" << std::endl;
    if( verbose )
      {
      std::cout << "Loading: " <<  inputs[0] << std::endl;
      }

    reader->SetFileName(inputs[0]);
    reader->Update();
    logfilt->SetInput(reader->GetOutput() );
    logfilt->Update();

    dup->SetInputImage(logfilt->GetOutput() );
    dup->Update();
    LogTensorImageType::Pointer average = dup->GetOutput();
    for( int i = 1; i < numberofinputs; ++i )
      {
      if( verbose )
        {
        std::cout << "Loading: " <<  inputs[i] << std::endl;
        }
      reader->SetFileName(inputs[i].c_str() );

      adder->SetInput1(average);
      adder->SetInput2(logfilt->GetOutput() );
      adder->Update();

      dup->SetInputImage(adder->GetOutput() );
      dup->Update();
      average = dup->GetOutput();

      }

    DivideImageFilter::Pointer divide = DivideImageFilter::New();
    divide->SetFunctor(PixelDivider<LogTensorPixelType>(numberofinputs) );
    divide->SetInput(average);
    divide->Update();

    typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
    ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
    expf->SetInput(divide->GetOutput() );

    if( !doubleDTI )
      {
      typedef itk::DiffusionTensor3D<float> TensorFloatPixelType;
      typedef itk::Image<TensorFloatPixelType, 3> TensorFloatImageType;
      typedef itk::CastImageFilter< TensorImageType, TensorFloatImageType > CastDTIFilterType ;
      CastDTIFilterType::Pointer castFilter = CastDTIFilterType::New() ;
      castFilter->SetInput( expf->GetOutput() ) ;
      typedef itk::ImageFileWriter<TensorFloatImageType> TensorFileWriterType;
      TensorFileWriterType::Pointer tensorWriter = TensorFileWriterType::New();
      tensorWriter->SetFileName(tensorOutput.c_str());
      tensorWriter->SetInput(castFilter->GetOutput());
      tensorWriter->SetUseCompression(true);
      tensorWriter->Update();
      }
    else
      {
      typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;
      TensorFileWriterType::Pointer twrit = TensorFileWriterType::New();
      twrit->SetUseCompression(true);
      twrit->SetInput(expf->GetOutput());
      twrit->SetFileName(tensorOutput);
      twrit->Update();
      }
    }
  else
    {
    std::cout << "At least one tensor field has to be specified" << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
