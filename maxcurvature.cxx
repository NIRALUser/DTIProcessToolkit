#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkNthElementImageAdaptor.h>
#include <itkAbsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkCastImageFilter.h>

#include <vnl/vnl_vector_fixed.h>


int main(int argc, const char* argv[])
{
  if(argc != 4) 
    {
    std::cerr << "Usage: " << argv[0] << " <infile> <outfile> <scale>" << std::endl;
    std::cerr << "Produces the image of the maximum curvature." << std::endl;
    return -1;
    
    }

  typedef unsigned short PixelType;
  typedef double FloatPixelType;
  const int DIM = 3;
  typedef itk::SymmetricSecondRankTensor<double, DIM> HessianPixelType;
  typedef itk::Vector<double, DIM> EigenValuePixelType;

  typedef itk::Image<PixelType, DIM> ImageType;
  typedef itk::Image<FloatPixelType, DIM> FloatImageType;

  typedef itk::Image<HessianPixelType, DIM> HessianImageType;
  typedef itk::Image<EigenValuePixelType, DIM> EigenValueImageType;

  typedef itk::ImageFileReader<ImageType> FileReaderType;
 
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(argv[1]);

  typedef itk::HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;
  HessianFilterType::Pointer hessian = HessianFilterType::New();
  hessian->SetInput(reader->GetOutput());
  hessian->SetSigma(atof(argv[3]));

  typedef itk::SymmetricEigenAnalysisImageFilter<HessianImageType,EigenValueImageType> EigenAnalysisFilterType;
  EigenAnalysisFilterType::Pointer eigen = EigenAnalysisFilterType::New();
  eigen->SetInput(hessian->GetOutput());
  eigen->OrderEigenValuesBy(EigenAnalysisFilterType::FunctorType::OrderByValue);
  eigen->SetDimension(3);

  eigen->Update();
  typedef itk::NthElementImageAdaptor<EigenValueImageType,FloatPixelType> ElementSelectAdaptorType;
  ElementSelectAdaptorType::Pointer elementSelect = ElementSelectAdaptorType::New();
  elementSelect->SetImage(eigen->GetOutput());
  elementSelect->SelectNthElement(0);

//   typedef itk::AbsImageFilter<ElementSelectAdaptorType,FloatImageType> AbsImageType;
//   AbsImageType::Pointer absfilt = AbsImageType::New();
//   absfilt->SetInput(elementSelect);

  typedef itk::CastImageFilter<ElementSelectAdaptorType,FloatImageType> CastFilter2Type;
  CastFilter2Type::Pointer cast2 = CastFilter2Type::New();
  cast2->SetInput(elementSelect);

  typedef itk::ShiftScaleImageFilter<FloatImageType,FloatImageType> ScaleImageType;
  ScaleImageType::Pointer scale = ScaleImageType::New();
  scale->SetScale(-10.0);
  scale->SetShift(0.0);
  scale->SetInput(cast2->GetOutput());

  scale->Update();
  std::cout << "scale: " << scale->GetOutput()->GetSpacing() << std::endl;

  typedef itk::ThresholdImageFilter<FloatImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  threshold->SetInput(scale->GetOutput());
  threshold->ThresholdBelow(0);

  typedef itk::CastImageFilter<FloatImageType,ImageType> CastFilterType;
  CastFilterType::Pointer cast = CastFilterType::New();
  cast->SetInput(threshold->GetOutput());

  typedef itk::ImageFileWriter<ImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetInput(cast->GetOutput());
  writer->SetFileName(argv[2]);

  try
    {
    eigen->Update();
    writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e << std::endl;
    return -1;
    }


  return 0;
}
