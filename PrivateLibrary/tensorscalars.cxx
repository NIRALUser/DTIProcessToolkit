#include "tensorscalars.h"

#include <itkMetaDataObject.h>

// Filters
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// My ITK Filters
#include "itkVectorMaskNegatedImageFilter.h"
#include "itkTensorMeanDiffusivityImageFilter.h"
#include "itkTensorColorFAImageFilter.h"
#include "itkTensorNegativeEigenValueImageFilter.h"
#include "itkTensorPrincipalEigenvectorImageFilter.h"
#include "itkVectorClosestDotProductImageFilter.h"
#include "itkTensorFAGradientImageFilter.h"
#include "itkTensorRotateImageFilter.h"

// Global constants
const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

void createFA(TensorImageType::Pointer timg,      // Tensor image
              const std::string &filename,   // Output filename
              bool scale)                    // if true write as int,
                                             // if false write as float
{
  typedef itk::TensorFractionalAnisotropyImageFilter<TensorImageType,RealImageType> FAFilterType;
  FAFilterType::Pointer fafilter = FAFilterType::New();
  fafilter->SetInput(timg);
    
  // If scale option set scale fa, and write out an integer image
  if(scale)
    {
//         unsigned int scale = vm["fa-scale"].as<unsigned int>();

    typedef itk::ShiftScaleImageFilter<RealImageType,IntImageType> ShiftScaleFilterType;
    ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
    scalefilter->SetInput(fafilter->GetOutput());
    scalefilter->SetShift(0);
    scalefilter->SetScale(10000);

    typedef itk::ImageFileWriter<IntImageType> IntImageWriterType;
    IntImageWriterType::Pointer intwriter = IntImageWriterType::New();
    intwriter->SetUseCompression(true);
    intwriter->SetInput(scalefilter->GetOutput());
    intwriter->SetFileName(filename.c_str());
    intwriter->Update();
    }
  // otherwise write out fa as float image w/ pixel \in [0,1]
  else
    {
    typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
    RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
    realwriter->SetUseCompression(true);
    realwriter->SetInput(fafilter->GetOutput());
    realwriter->SetFileName(filename.c_str());
    realwriter->Update();
         
    }

}

void createFAGradient(TensorImageType::Pointer timg, // Tensor image
                      const std::string &filename,   // Output filename
                      double sigma)                  // sigma
{
  typedef itk::TensorFAGradientImageFilter<double> FAGradientImageFilter;

  FAGradientImageFilter::Pointer fagradfilter = FAGradientImageFilter::New();
  fagradfilter->SetInput(timg);
  fagradfilter->SetSigma(sigma);

  typedef itk::ImageFileWriter<GradientImageType> GradientImageWriter;
  GradientImageWriter::Pointer gradwriter = GradientImageWriter::New();
  gradwriter->SetUseCompression(true);
  gradwriter->SetInput(fagradfilter->GetOutput());
  gradwriter->SetFileName(filename.c_str());
  gradwriter->Update();
}

void createFAGradMag(TensorImageType::Pointer timg,      // Tensor image
                     const std::string &filename,        // Output
                                                         // filename
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

  typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
  RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
  realwriter->SetUseCompression(true);
  realwriter->SetInput(gradmag->GetOutput());
  realwriter->SetFileName(filename.c_str());
  realwriter->Update();
  

}


void createColorFA(TensorImageType::Pointer timg,      // Tensor image
              const std::string &filename)   // Output filename
{
  typedef itk::RGBPixel<unsigned char> RGBPixel;
  typedef itk::Image<RGBPixel,3> RGBImageType;
  typedef itk::TensorColorFAImageFilter<TensorImageType,RGBImageType> FAFilterType;
  FAFilterType::Pointer fafilter = FAFilterType::New();
  fafilter->SetInput(timg);
    
  typedef itk::ImageFileWriter<RGBImageType> IntImageWriterType;
  IntImageWriterType::Pointer intwriter = IntImageWriterType::New();
  intwriter->SetUseCompression(true);
  intwriter->SetInput(fafilter->GetOutput());
  intwriter->SetFileName(filename.c_str());
  intwriter->Update();
}

void createPrincipalEigenvector(TensorImageType::Pointer timg,      // Tensor image
                                const std::string &eigvecname,    // Output
                                                                  // principal
                                                                  // eigenvector name
                                const std::string &dotprodname,   // Ouput
                                                                  // closest
                                                                  // dotproduct
                                                                  // name
                                GradientListType::Pointer gradientlist)
                                
{
    itk::MetaDataDictionary & dict = timg->GetMetaDataDictionary();

    if(dict.HasKey(NRRD_MEASUREMENT_KEY))
      {
      // measurement frame
      vnl_matrix<double> mf(3,3);
      // imaging frame
      vnl_matrix<double> imgf(3,3);
      
      std::vector<std::vector<double> > nrrdmf;
      itk::ExposeMetaData<std::vector<std::vector<double> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);
      
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
      
      itk::EncapsulateMetaData<std::vector<std::vector<double> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);
      
      typedef itk::TensorRotateImageFilter<TensorImageType, TensorImageType, double> TensorRotateFilterType;
      TensorRotateFilterType::Pointer trotate = TensorRotateFilterType::New();
      trotate->SetInput(timg);
      trotate->SetRotation(vnl_svd<double>(imgf).inverse()*mf);
      trotate->Update();
      timg = trotate->GetOutput();
      
    }

  typedef itk::CovariantVector<double, 3> VectorPixel;
  typedef itk::Image<VectorPixel, 3> VectorImageType;
  typedef itk::TensorPrincipalEigenvectorImageFilter<TensorImageType,VectorImageType> PrincipalEigenvectorFilterType;

  PrincipalEigenvectorFilterType::Pointer ppdfilter = PrincipalEigenvectorFilterType::New();
  ppdfilter->SetInput(timg);

  if(eigvecname != "")
    {
    typedef itk::ImageFileWriter<VectorImageType> VectorImageWriterType;
    VectorImageWriterType::Pointer vectorwriter = VectorImageWriterType::New();
    vectorwriter->SetInput(ppdfilter->GetOutput());
    vectorwriter->SetFileName(eigvecname.c_str());
    try
      {
      vectorwriter->Update();
      }
    catch (itk::ExceptionObject & e)
      {
      std::cerr << "Couldn't write principal eigenvector field" << std::endl;
      std::cerr << e << std::endl;
      }
    }

  if(dotprodname != "")
    {
 
    typedef itk::VectorClosestDotProductImageFilter<VectorImageType, RealImageType> ClosestDotProductFilterType;
    ClosestDotProductFilterType::Pointer cdpfilter = ClosestDotProductFilterType::New();
    cdpfilter->SetInput(ppdfilter->GetOutput());
    cdpfilter->SetGradientList(gradientlist);

    typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
    RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
    realwriter->SetUseCompression(true);
    realwriter->SetInput(cdpfilter->GetOutput());
    realwriter->SetFileName(dotprodname.c_str());

    try
      {
      realwriter->Update();
      }
    catch (itk::ExceptionObject & e)
      {
      std::cerr << "Couldn't write closest dot product image" << std::endl;
      std::cerr << e << std::endl;
      }
    } // end if(dotprodname != "")
}

void createMD(TensorImageType::Pointer timg,
              const std::string &filename,
              bool scale)
{
  typedef itk::TensorMeanDiffusivityImageFilter<TensorImageType,RealImageType> MDFilterType;
  MDFilterType::Pointer mdfilter = MDFilterType::New();
  mdfilter->SetInput(timg);
    
  // If scale option set scale md, and write out an integer image
  if(scale)
    {
    typedef itk::ShiftScaleImageFilter<RealImageType,IntImageType> ShiftScaleFilterType;
    ShiftScaleFilterType::Pointer scalefilter = ShiftScaleFilterType::New();
    scalefilter->SetInput(mdfilter->GetOutput());
    scalefilter->SetShift(0);
    scalefilter->SetScale(100000);

    typedef itk::ImageFileWriter<IntImageType> IntImageWriterType;
    IntImageWriterType::Pointer intwriter = IntImageWriterType::New();
    intwriter->SetUseCompression(true);
    intwriter->SetInput(scalefilter->GetOutput());
    intwriter->SetFileName(filename.c_str());
    intwriter->Update();
    }
  // otherwise write out md as float image w/ pixel \in [0,1]
  else
    {
    typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
    RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
    realwriter->SetUseCompression(true);
    realwriter->SetInput(mdfilter->GetOutput());
    realwriter->SetFileName(filename.c_str());
    realwriter->Update();
         
    }
}
void createNegativeEigenValueLabel(TensorImageType::Pointer timg,
                                   const std::string &filename)
{
  typedef itk::TensorNegativeEigenValueImageFilter<TensorImageType,
    LabelImageType> NegEigLabelFilterType;

  NegEigLabelFilterType::Pointer negeigdetect = NegEigLabelFilterType::New();
  negeigdetect->SetInput(timg);

  
  // TODO: wrap in try catch block
  typedef itk::ImageFileWriter<LabelImageType> LabelImageWriterType;
  LabelImageWriterType::Pointer labelwriter = LabelImageWriterType::New();
  labelwriter->SetUseCompression(true);
  labelwriter->SetInput(negeigdetect->GetOutput());
  labelwriter->SetFileName(filename.c_str());
  labelwriter->Update();
  
}

