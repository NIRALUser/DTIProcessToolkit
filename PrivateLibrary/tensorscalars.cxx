#include "tensorscalars.h"

#include <itkMetaDataObject.h>

// Filters
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkFastSymmetricEigenAnalysisImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

// My ITK Filters
#include "itkVectorMaskNegatedImageFilter.h"
#include "itkTensorMeanDiffusivityImageFilter.h"
#include "itkTensorFrobeniusNormImageFilter.h"
#include "itkTensorColorFAImageFilter.h"
#include "itkTensorNegativeEigenValueImageFilter.h"
#include "itkTensorPrincipalEigenvectorImageFilter.h"
#include "itkVectorClosestDotProductImageFilter.h"
#include "itkTensorFAGradientImageFilter.h"
#include "itkTensorRotateImageFilter.h"

// Global constants
const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

template<>
itk::Image<double, 3>::Pointer createFA<double>(TensorImageType::Pointer timg) // Tensor image
{
  typedef itk::TensorFractionalAnisotropyImageFilter<TensorImageType,RealImageType> FAFilterType;
  FAFilterType::Pointer fafilter = FAFilterType::New();
  fafilter->SetInput(timg);
  fafilter->Update();
    
  return fafilter->GetOutput();
}

template<>
itk::Image<unsigned short, 3>::Pointer createFA<unsigned short>(TensorImageType::Pointer timg)      // Tensor image
{
  RealImageType::Pointer realfa = createFA<double>(timg);
    
  typedef itk::ShiftScaleImageFilter<RealImageType,IntImageType> ShiftScaleFilterType;
  ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
  scalefilter->SetInput(realfa);
  scalefilter->SetShift(0);
  scalefilter->SetScale(10000);
  scalefilter->Update();

  return scalefilter->GetOutput();
}

template<>
itk::Image<double, 3>::Pointer createMD<double>(TensorImageType::Pointer timg) // Tensor image
{
  typedef itk::TensorMeanDiffusivityImageFilter<TensorImageType,RealImageType> MDFilterType;
  MDFilterType::Pointer mdfilter = MDFilterType::New();
  mdfilter->SetInput(timg);
  mdfilter->Update();
    
  return mdfilter->GetOutput();
}

template<>
itk::Image<unsigned short, 3>::Pointer createMD<unsigned short>(TensorImageType::Pointer timg)      // Tensor image
{
  RealImageType::Pointer realmd = createMD<double>(timg);
    
  typedef itk::ShiftScaleImageFilter<RealImageType,IntImageType> ShiftScaleFilterType;
  ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
  scalefilter->SetInput(realmd);
  scalefilter->SetShift(0);
  scalefilter->SetScale(100000);
  scalefilter->Update();

  return scalefilter->GetOutput();
}

template<>
itk::Image<double, 3>::Pointer createLambda<double>(TensorImageType::Pointer timg, // Tensor image
                                                    EigenValueIndex lambdaind) // Lambda index
{
  // Not really a deformation image output jsut a 3-vector of doubles.
  typedef itk::FastSymmetricEigenAnalysisImageFilter<TensorImageType,DeformationImageType> LambdaFilterType;
  LambdaFilterType::Pointer lambdafilter = LambdaFilterType::New();
  lambdafilter->SetInput(timg);
  lambdafilter->OrderEigenValuesBy(LambdaFilterType::FunctorType::OrderByValue);
  lambdafilter->Update();
    
  typedef itk::VectorIndexSelectionCastImageFilter<LambdaFilterType::OutputImageType, RealImageType> ElementSelectAdaptorType;
  ElementSelectAdaptorType::Pointer elementSelect = ElementSelectAdaptorType::New();
  elementSelect->SetInput(lambdafilter->GetOutput());
  // Reverse semanatics of lambda_1 from ITK.  
  // In our convention lambda_1 is the largest eigenvalue whereas in
  // ITK its the smallest
  elementSelect->SetIndex(2 - lambdaind); 
  elementSelect->Update();

  return elementSelect->GetOutput();
}

template<>
itk::Image<unsigned short, 3>::Pointer createLambda<unsigned short>(TensorImageType::Pointer timg,      // Tensor image
                                                                    EigenValueIndex lambdaind) // Lambda index
{
  RealImageType::Pointer reallambda = createLambda<double>(timg, lambdaind);
    
  typedef itk::ShiftScaleImageFilter<RealImageType,IntImageType> ShiftScaleFilterType;
  ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
  scalefilter->SetInput(reallambda);
  scalefilter->SetShift(0);
  scalefilter->SetScale(100000);
  scalefilter->Update();

  return scalefilter->GetOutput();
}


template<>
itk::Image<double, 3>::Pointer createFro<double>(TensorImageType::Pointer timg) // Tensor image
{
  typedef itk::TensorFrobeniusNormImageFilter<TensorImageType,RealImageType> MDFilterType;
  MDFilterType::Pointer mdfilter = MDFilterType::New();
  mdfilter->SetInput(timg);
  mdfilter->Update();
    
  return mdfilter->GetOutput();
}

template<>
itk::Image<unsigned short, 3>::Pointer createFro<unsigned short>(TensorImageType::Pointer timg)      // Tensor image
{
  RealImageType::Pointer realfro = createFro<double>(timg);
    
  typedef itk::ShiftScaleImageFilter<RealImageType,IntImageType> ShiftScaleFilterType;
  ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
  scalefilter->SetInput(realfro);
  scalefilter->SetShift(0);
  scalefilter->SetScale(100000);
  scalefilter->Update();

  return scalefilter->GetOutput();
}

GradientImageType::Pointer createFAGradient(TensorImageType::Pointer timg, // Tensor image
                                            double sigma)                  // sigma
{
  typedef itk::TensorFAGradientImageFilter<double> FAGradientImageFilter;

  FAGradientImageFilter::Pointer fagradfilter = FAGradientImageFilter::New();
  fagradfilter->SetInput(timg);
  fagradfilter->SetSigma(sigma);

  return fagradfilter->GetOutput();
}

RealImageType::Pointer createFAGradMag(TensorImageType::Pointer timg,      // Tensor image
                                       double sigma)
{
  typedef itk::TensorFractionalAnisotropyImageFilter<TensorImageType,RealImageType> FAFilterType;
  FAFilterType::Pointer fafilter = FAFilterType::New();
  fafilter->SetInput(timg);
    
  // If scale option set scale fa, and write out an integer image


  typedef itk::ShiftScaleImageFilter<RealImageType,RealImageType> ShiftScaleFilterType;
  ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
  scalefilter->SetInput(fafilter->GetOutput());
  scalefilter->SetShift(0);
  scalefilter->SetScale(10000);

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<RealImageType,RealImageType> GradMagFilterType;
  GradMagFilterType::Pointer gradmag = GradMagFilterType::New();
  gradmag->SetInput(scalefilter->GetOutput());
  gradmag->SetSigma(sigma);
  gradmag->SetNormalizeAcrossScale(false);
  gradmag->Update();

  return gradmag->GetOutput();
}


RGBImageType::Pointer createColorFA(TensorImageType::Pointer timg)      // Tensor image
{
  typedef itk::RGBPixel<unsigned char> RGBPixel;
  typedef itk::Image<RGBPixel,3> RGBImageType;
  typedef itk::TensorColorFAImageFilter<TensorImageType,RGBImageType> FAFilterType;
  FAFilterType::Pointer fafilter = FAFilterType::New();
  fafilter->SetInput(timg);
  fafilter->Update();
  
  return fafilter->GetOutput();
}

GradientImageType::Pointer createPrincipalEigenvector(TensorImageType::Pointer timg) // Tensor image
{
    itk::MetaDataDictionary & dict = timg->GetMetaDataDictionary();

    if(dict.HasKey(NRRD_MEASUREMENT_KEY))
      {
      // measurement frame
      vnl_matrix<double> mf(3,3);
      // imaging frame
      vnl_matrix<double> imgf(3,3);
      
      std::vector<std::vector<double> > nrrdmf;
      itk::ExposeMetaData<std::vector<std::vector<double> > >(dict, NRRD_MEASUREMENT_KEY, nrrdmf);
      
      imgf = timg->GetDirection().GetVnlMatrix();
      for(unsigned int i = 0; i < 3; ++i)
        {
        for(unsigned int j = 0; j < 3; ++j)
          {
          mf(i,j) = nrrdmf[i][j];
          
          if(i == j)
            nrrdmf[i][j] = 1.0;
          else
            nrrdmf[i][j] = 0.0;
          }
        }
      
      itk::EncapsulateMetaData<std::vector<std::vector<double> > >(dict, NRRD_MEASUREMENT_KEY, nrrdmf);
      
      typedef itk::TensorRotateImageFilter<TensorImageType, TensorImageType, double> TensorRotateFilterType;
      TensorRotateFilterType::Pointer trotate = TensorRotateFilterType::New();
      trotate->SetInput(timg);
      trotate->SetRotation(vnl_svd<double>(imgf).inverse() * mf);
      trotate->Update();
      timg = trotate->GetOutput();
      
    }

  typedef itk::CovariantVector<double, 3> VectorPixel;
  typedef itk::Image<VectorPixel, 3> VectorImageType;
  typedef itk::TensorPrincipalEigenvectorImageFilter<TensorImageType,VectorImageType> PrincipalEigenvectorFilterType;

  PrincipalEigenvectorFilterType::Pointer ppdfilter = PrincipalEigenvectorFilterType::New();
  ppdfilter->SetInput(timg);
  ppdfilter->Update();
  return ppdfilter->GetOutput();
}

LabelImageType::Pointer createNegativeEigenValueLabel(TensorImageType::Pointer timg)
{
  typedef itk::TensorNegativeEigenValueImageFilter<TensorImageType,
    LabelImageType> NegEigLabelFilterType;

  NegEigLabelFilterType::Pointer negeigdetect = NegEigLabelFilterType::New();
  negeigdetect->SetInput(timg);
  negeigdetect->Update();
 
  return negeigdetect->GetOutput();
}

