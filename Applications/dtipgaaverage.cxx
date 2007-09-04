#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_matrix_inverse.h>

#include <itkAffineTransform.h>
#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVector.h>

//#include "itkLinearInterpolateTensorImageImageFunction.h"
#include "itkNaryTensorAverageImageFilter.h"

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

  std::vector<TensorFileReader::Pointer> readers(sources.size());

  const int n = sources.size();

  typedef itk::NaryTensorAverageImageFilter<TensorImageType,TensorImageType> NaryAverageFilterType;
  NaryAverageFilterType::Pointer averagefilt = NaryAverageFilterType::New();
  averagefilt->SetNumberOfThreads(1);

  for(int i = 0; i < n; ++i)
    {
    readers[i] = TensorFileReader::New();
    readers[i]->SetFileName(sources[i].c_str());
    readers[i]->Update();
    
    averagefilt->SetInput(i,readers[i]->GetOutput());
    }
  averagefilt->Update();


  std::cout << "Writing result" << std::endl;
  typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;
  TensorFileWriterType::Pointer twrit = TensorFileWriterType::New();
  twrit->SetInput(averagefilt->GetOutput());
  twrit->SetFileName(ofile.c_str());
  twrit->Update();

}

