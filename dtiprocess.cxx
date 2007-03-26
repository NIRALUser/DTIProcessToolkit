// STL includes
#include <string>
#include <iostream>
#include <fstream>

// boost includes
#include <boost/program_options.hpp>

// ITK includes
// datastructures
#include <itkImage.h>
#include <itkVector.h>
#include <itkMetaDataObject.h>
#include <itkVector.h>
#include <itkCovariantVector.h>
#include <itkDiffusionTensor3D.h>
#include <itkAffineTransform.h>

// IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Filters
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkVectorNearestNeighborInterpolateImageFunction.h>

// VNL Includes
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>

// My ITK Filters
#include "itkVectorMaskNegatedImageFilter.h"
#include "itkTensorRotateImageFilter.h"
#include "itkTensorRotateFromDeformationFieldImageFilter.h"
#include "itkTensorRotateFromDeformationFieldPPDImageFilter.h"
#include "itkDeformationFieldJacobianFilter.h"
#include "itkLogEuclideanTensorImageFilter.h"
#include "itkExpEuclideanTensorImageFilter.h"
#include "itkTensorMeanDiffusivityImageFilter.h"
#include "itkTensorColorFAImageFilter.h"
#include "itkHFieldToDeformationFieldImageFilter.h"
#include "itkTensorNegativeEigenValueImageFilter.h"
#include "itkTensorPrincipalEigenvectorImageFilter.h"
#include "itkVectorClosestDotProductImageFilter.h"
#include "itkVectorBSplineInterpolateImageFunction.h"
#include "itkTensorFAGradientImageFilter.h"

// dtiprocess headers
#include "transforms.h"

const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

namespace po = boost::program_options;

// Define necessary types for images
typedef double RealType;
typedef double TransformRealType;
typedef unsigned char LabelType;
const unsigned int DIM = 3;

typedef unsigned short ScalarPixelType;
typedef itk::DiffusionTensor3D<double> TensorPixelType;
typedef itk::Vector<double,3> DeformationPixelType;
typedef itk::CovariantVector<double,3> GradientPixelTyep;

typedef itk::Image<TensorPixelType, DIM> TensorImageType;

typedef itk::Image<DeformationPixelType, DIM> DeformationImageType;
typedef itk::Image<GradientPixelTyep, DIM> GradientImageType;

typedef itk::Image<RealType, DIM> RealImageType;
typedef itk::Image<ScalarPixelType, DIM> IntImageType;
typedef itk::Image<LabelType, DIM> LabelImageType;

typedef TensorImageType::SizeType ImageSizeType;
typedef TensorImageType::SpacingType ImageSpacingType;

typedef itk::AffineTransform<TransformRealType,3> AffineTransformType;

typedef vnl_vector_fixed<double, 3> GradientType;
typedef itk::VectorContainer<unsigned int, GradientType> GradientListType;

enum InterpolationType {NearestNeighbor, Linear, Cubic};
enum TensorReorientationType {FiniteStrain, PreservationPrincipalDirection};

// derived output functions
void createFA(TensorImageType::Pointer, const std::string &, bool);
void createFAGradient(TensorImageType::Pointer, const std::string &,double);
void createFAGradMag(TensorImageType::Pointer, const std::string &,double);
void createColorFA(TensorImageType::Pointer, const std::string &);
void createPrincipalEigenvector(TensorImageType::Pointer, const std::string &, const std::string &, GradientListType::Pointer);
void createMD(TensorImageType::Pointer, const std::string &, bool);
void createNegativeEigenValueLabel(TensorImageType::Pointer, const std::string &);

// warping functions
void createROT(TensorImageType::Pointer, const std::string &, const std::string &);
void createWarp(TensorImageType::Pointer, const std::string &, const std::string &, const std::string &, TensorReorientationType, InterpolationType);

// Bad global variables.  TODO: remove these
bool VERBOSE=false;


void validate(boost::any& v,
              const std::vector<std::string>& values,
              InterpolationType* target_type,
              int)
{
  using namespace po;
  using boost::any;

  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s ==  "nearestneighbor")
    {
    v = any(NearestNeighbor);
    }
  else if (s == "linear")
    {
    v = any(Linear);
    } 
  else if (s == "cubic")
    {
    v = any(Cubic);
    }
  else
    {
    throw validation_error("Interpolation type invalid.  Only \"nearestneighbor\", \"linear\", and \"cubic\" allowed.");
    }
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              TensorReorientationType* target_type,
              int)
{
  using namespace po;
  using boost::any;
  // Make sure no previous assignment to 'a' was made.
  validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = validators::get_single_string(values);

  if(s == "fs")
    {
    v = any(FiniteStrain);
    }
  else if(s == "ppd")
    {
    v = any(PreservationPrincipalDirection);
    }
  else
    {
    throw validation_error("Reorientation type invalid.  Only \"fs\" or \"ppd\"");
    }
}


int main(int argc, char* argv[])
{

  // Read program options/configuration
  po::options_description config("Usage: dtiprocess tensor-image [options]");
  config.add_options()
    // General options
    ("help,h", "produce this help message")
    ("verbose,v","produces verbose output")

    // Derived outputs
    ("fa-output,f", po::value<std::string>(),"FA output file")
    ("fa-gradient-output", po::value<std::string>(),"FA gradient")
    ("sigma,s", po::value<double>()->default_value(2.0), "Scale of gradients")
    ("fa-gradmag-output", po::value<std::string>(),"FA gradient magnitude")
    ("color-fa-output,c", po::value<std::string>(),"Color FA output file")
    ("principal-eigenvector-output,V", po::value<std::string>(),"Principal Eigenvector of tensor field")
    ("closest-dotproduct-output,D", po::value<std::string>(), "Closes dot product of principal eigenvector to all gradient directions")
    ("md-output,m", po::value<std::string>(),"MD output file")
    ("negative-eig-output,n", po::value<std::string>(),"Negative eigenvalue detected output file")

    // derived output options
    ("scalar-float","Write scalar [FA,MD] as unscaled float.  Also causes FA to be unscaled [0..1].")
    ("double", "Writes output in double precisision")

    // tensor transformations
    // affine
    ("rot-output,r", po::value<std::string>(),"Rotated tensor output file.  Must also specify the dof file")
    ("dof-file,d", po::value<std::string>(), "DOF for rotation")

    // deformation
    ("deformation-output,w", po::value<std::string>(), "Warped tensor field based on a deformation field.  Must input h as \"h-field\" of transform")
    ("h-field,H", po::value<std::string>(), "HField for warp")
    ("inv-h-field,I", po::value<std::string>(), "Inverse HField for warp")

    // transformation options
    ("interpolation,i", po::value<InterpolationType>()->default_value(Linear,"linear"), "Interpolation type")
    ("reorientation", po::value<TensorReorientationType>()->default_value(FiniteStrain,"finite strain"), "Reorientation type")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("dti-image", po::value<std::string>(), "DTI tensor volume")
  ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("dti-image",1);

  po::variables_map vm;

  try
    {
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);     
    } 
  catch (const po::error &e)
    {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
    }

  // End option reading configuration

  // Display help if asked or program improperly called
  if(vm.count("help") || !vm.count("dti-image"))
    {
    std::cout << config << std::endl;
    if(vm.count("help"))
      return EXIT_SUCCESS;
    else
      return EXIT_FAILURE;
    }

  bool scale = true;
  if(vm.count("scalar-float"))
    {
    scale = false;
    }

  if(vm.count("verbose"))
    {
    VERBOSE = true;
    }

  typedef itk::ImageFileReader<TensorImageType> FileReaderType;
  FileReaderType::Pointer dtireader = FileReaderType::New();
  dtireader->SetFileName(vm["dti-image"].as<std::string>().c_str());

  try
    {
    dtireader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << e <<std::endl;
    return EXIT_FAILURE;
    }

  bool readb0 = false;
  double b0 = 0;

  GradientListType::Pointer gradientContainer = GradientListType::New();

  itk::MetaDataDictionary & dict = dtireader->GetOutput()->GetMetaDataDictionary();
  std::vector<std::string> keys = dict.GetKeys();
  for(std::vector<std::string>::const_iterator it = keys.begin();
      it != keys.end(); ++it)
    {
    std::string value;
    if( it->find("DWMRI_b-value") != std::string::npos)
      {
      std::string t;
      itk::ExposeMetaData<std::string>(dict, *it, t);
      readb0 = true;
      b0 = atof(t.c_str());
      }
    else if( it->find("DWMRI_gradient") != std::string::npos)
      {
      std::string value;

      itk::ExposeMetaData<std::string>(dict, *it, value);
      std::istringstream iss(value);
      GradientType g;
      iss >> g[0] >> g[1] >> g[2];

      if(g.two_norm() > 0.0)
        g = g / g.two_norm();

      unsigned int ind;
      std::string temp = it->substr(it->find_last_of('_')+1);
      ind = atoi(temp.c_str());
      
      gradientContainer->InsertElement(ind,g);
      }
    else if( it->find("DWMRI_NEX") != std::string::npos)
      {
      std::string numrepstr;

      itk::ExposeMetaData<std::string>(dict, *it, numrepstr);
      unsigned int numreps = atoi(value.c_str());

      std::string indtorepstr = it->substr(it->find_last_of('_')+1);
      unsigned int indtorep =  atoi(indtorepstr.c_str());

      GradientType g = gradientContainer->GetElement(indtorep);

      for(unsigned int i = indtorep+1; i < indtorep+numreps-1; i++)
        gradientContainer->InsertElement(i,g);
      }

    }

  // Debugging
  if(VERBOSE)
    {
    std::cout << "Interpolation type: " << vm["interpolation"].as<InterpolationType>() << std::endl;
    std::cout << "reorientation type: " << vm["reorientation"].as<TensorReorientationType>() << std::endl;
    }


  double sigma = vm["sigma"].as<double>();

  // Compute FA image
  if(vm.count("fa-output"))
    {
    createFA(dtireader->GetOutput(),vm["fa-output"].as<std::string>(),scale);
    }

  if(vm.count("fa-gradient-output"))
    {
    createFAGradient(dtireader->GetOutput(),vm["fa-gradient-output"].as<std::string>(),sigma);
    }

  if(vm.count("fa-gradmag-output"))
    {
    createFAGradMag(dtireader->GetOutput(),vm["fa-gradmag-output"].as<std::string>(),sigma);
    }

  if(vm.count("color-fa-output"))
    {
    createColorFA(dtireader->GetOutput(),vm["color-fa-output"].as<std::string>());
    }

  if(vm.count("principal-eigenvector-output") || vm.count("closest-dotproduct-output"))
    {
    const std::string peo(vm.count("principal-eigenvector-output") ?
                          vm["principal-eigenvector-output"].as<std::string>() :
                          "");

    const std::string cdo(vm.count("closest-dotproduct-output") ?
                          vm["closest-dotproduct-output"].as<std::string>() :
                          "");

    createPrincipalEigenvector(dtireader->GetOutput(),peo,cdo,gradientContainer);
    }

  if(vm.count("md-output"))
    {
    createMD(dtireader->GetOutput(),vm["md-output"].as<std::string>(),scale);
    }

  if(vm.count("negative-eig-output"))
    {
    createNegativeEigenValueLabel(dtireader->GetOutput(),vm["negative-eig-output"].as<std::string>());
    }


  if(vm.count("rot-output"))
    {
    if(!vm.count("dof-file"))
      {
      std::cerr << "Tensor rotation requested, but dof file not specified" << std::endl;
      return EXIT_FAILURE;
      }
    createROT(dtireader->GetOutput(),
              vm["rot-output"].as<std::string>(),
              vm["dof-file"].as<std::string>());
    }
   

  if(vm.count("deformation-output"))
    {
    if(!vm.count("h-field") || !vm.count("inv-h-field"))
      {
      std::cerr << "Deformation field info not fully specified" << std::endl;
      return EXIT_FAILURE;
      }
    createWarp(dtireader->GetOutput(),
               vm["deformation-output"].as<std::string>(),
               vm["h-field"].as<std::string>(),
               vm["inv-h-field"].as<std::string>(),
               vm["reorientation"].as<TensorReorientationType>(),
               vm["interpolation"].as<InterpolationType>());
    }

  return EXIT_SUCCESS;
}


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
    intwriter->SetInput(scalefilter->GetOutput());
    intwriter->SetFileName(filename.c_str());
    intwriter->Update();
    }
  // otherwise write out fa as float image w/ pixel \in [0,1]
  else
    {
    typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
    RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
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
    intwriter->SetInput(scalefilter->GetOutput());
    intwriter->SetFileName(filename.c_str());
    intwriter->Update();
    }
  // otherwise write out md as float image w/ pixel \in [0,1]
  else
    {
    typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
    RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
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
  labelwriter->SetInput(negeigdetect->GetOutput());
  labelwriter->SetFileName(filename.c_str());
  labelwriter->Update();
  
  }

void createROT(TensorImageType::Pointer timg, 
               const std::string &ofile,
               const std::string &doffile)
{
  RViewTransform<TransformRealType> dof(readDOFFile<TransformRealType>(doffile));
  AffineTransformType::Pointer transform = 
    createITKAffine(dof,
                    timg->GetLargestPossibleRegion().GetSize(),
                    timg->GetSpacing());
  vnl_matrix<TransformRealType> R = 
    getInverseRotation(transform);

  if(VERBOSE)
    std::cout << "Rotating" << std::endl;
  typedef itk::TensorRotateImageFilter<TensorImageType,TensorImageType,TransformRealType> TensorRotateType;
  TensorRotateType::Pointer trotfilt = TensorRotateType::New();
  trotfilt->SetRotation(R);
  trotfilt->SetInput(timg);
  trotfilt->Update();

  if(VERBOSE)
    std::cout << "Taking Log" << std::endl;
  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType LogTensorImageType;
  LogEuclideanFilter::Pointer logf = LogEuclideanFilter::New();
  logf->SetInput(trotfilt->GetOutput()); 
  logf->Update();

  if(VERBOSE)
    std::cout << "Resampling" << std::endl;
  typedef itk::VectorResampleImageFilter<
    LogTensorImageType,LogTensorImageType> ResampleFilterType; // TODO: Should accept precision
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  // Set interpolator
  typedef itk::VectorBSplineInterpolateImageFunction<LogTensorImageType,double,double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  resampler->SetInterpolator(interpolator);

  resampler->SetInput(logf->GetOutput());
  resampler->SetTransform(transform);

  LogTensorImageType::Pointer logim = logf->GetOutput();
  resampler->SetSize( logim->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( logim->GetOrigin() );
  resampler->SetOutputSpacing( logim->GetSpacing() );
  
  typedef LogTensorImageType::PixelType LogPixelType;
  LogPixelType def(0.0);
  def[0] = -10;
  def[3] = -10;
  def[5] = -10;

  resampler->SetDefaultPixelValue(def);

  resampler->Update();

  if(VERBOSE)
     std::cout << "Taking Exp" << std::endl;
  typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
  ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
  expf->SetInput(resampler->GetOutput());
  expf->Update();

  if(VERBOSE)
    std::cout << "Writing result" << std::endl;
  typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;
  TensorFileWriterType::Pointer twrit = TensorFileWriterType::New();
  twrit->SetInput(expf->GetOutput());
  twrit->SetFileName(ofile.c_str());
  twrit->Update();
    
}

void createWarp(TensorImageType::Pointer timg,
                const std::string &ofile,
                const std::string &warpfile,
                const std::string &invwarpfile,
                TensorReorientationType reorientationtype,
                InterpolationType interpolationtype)
{
  // Read deformation field
  typedef itk::ImageFileReader<DeformationImageType> DeformationImageReader;
  DeformationImageReader::Pointer defreader = DeformationImageReader::New();
  defreader->SetFileName(warpfile.c_str());

  typedef itk::HFieldToDeformationFieldImageFilter<DeformationImageType> DeformationConvertType;
  DeformationConvertType::Pointer defconv = DeformationConvertType::New();
  defconv->SetInput(defreader->GetOutput());
//  defconv->SetSpacing(timg->GetSpacing());
  defconv->Update();

  // Read inverse deformation field
  DeformationImageReader::Pointer invreader = DeformationImageReader::New();
  invreader->SetFileName(invwarpfile.c_str());

  DeformationConvertType::Pointer invconv = DeformationConvertType::New();
  invconv->SetInput(invreader->GetOutput());

  // Compute jacobian of inverse deformation field
  typedef itk::DeformationFieldJacobianFilter<DeformationImageType,float> JacobianFilterType;
  typedef JacobianFilterType::OutputImageType JacobianImageType;
  JacobianFilterType::Pointer jacobian = JacobianFilterType::New();
  jacobian->SetInput(invconv->GetOutput());

  // Rotate tensor based on inverse deformation field
  typedef itk::InPlaceImageFilter<
    TensorImageType,
    TensorImageType> TensorRotateImageFilterBaseType;

  typedef itk::TensorRotateFromDeformationFieldImageFilter<
    TensorImageType,
    JacobianImageType,
    TensorImageType> TensorFSRotateImageFilterType;

  typedef itk::TensorRotateFromDeformationFieldPPDImageFilter<
    TensorImageType,
    JacobianImageType,
    TensorImageType> TensorPPDRotateImageFilterType;

  TensorRotateImageFilterBaseType::Pointer rotate;
  if(reorientationtype == FiniteStrain)
    {
    TensorFSRotateImageFilterType::Pointer fsrotate;
    fsrotate = TensorFSRotateImageFilterType::New();

    fsrotate->SetNumberOfThreads(1);
    fsrotate->SetInput1(timg);
    
    fsrotate->SetInput2(jacobian->GetOutput());

    rotate = fsrotate;
    }
  else if(reorientationtype == PreservationPrincipalDirection)
    {
    TensorPPDRotateImageFilterType::Pointer fsrotate;
    fsrotate = TensorPPDRotateImageFilterType::New();

    fsrotate->SetNumberOfThreads(1);
    fsrotate->SetInput1(timg);
    
    fsrotate->SetInput2(jacobian->GetOutput());
    rotate = fsrotate;
    }



  if(VERBOSE)
    std::cout << "Taking Log" << std::endl;
  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType LogTensorImageType;
  LogEuclideanFilter::Pointer logf = LogEuclideanFilter::New();
  logf->SetInput(rotate->GetOutput()); 
  logf->Update();

  typedef itk::WarpVectorImageFilter<LogTensorImageType,LogTensorImageType,DeformationImageType>
    WarpImageFilterType;
  WarpImageFilterType::Pointer warp = WarpImageFilterType::New();
  
  if(interpolationtype == Cubic)
    {
    typedef itk::VectorBSplineInterpolateImageFunction<LogTensorImageType,double,double> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    warp->SetInterpolator(interpolator);
    }
  else if (interpolationtype == Linear)
    {
    typedef itk::VectorLinearInterpolateImageFunction<LogTensorImageType,double> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    warp->SetInterpolator(interpolator);
    }  
  else if (interpolationtype == NearestNeighbor)
    {
    typedef itk::VectorNearestNeighborInterpolateImageFunction<LogTensorImageType,double> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    warp->SetInterpolator(interpolator);
    }

  warp->SetInput(logf->GetOutput());
  warp->SetDeformationField(defconv->GetOutput());
  warp->SetOutputSpacing(logf->GetOutput()->GetSpacing());
  warp->SetOutputOrigin(logf->GetOutput()->GetOrigin());
  warp->Update();

  typedef LogTensorImageType::PixelType LogPixelType;
  LogPixelType def(0.0);
  def[0] = -1e10;
  def[3] = -1e10;
  def[5] = -1e10;

  warp->SetEdgePaddingValue(def);
  warp->Update();

  if(VERBOSE)
     std::cout << "Taking Exp" << std::endl;
  typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
  ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
  expf->SetInput(warp->GetOutput());
  expf->Update();

  if(VERBOSE)
    std::cout << "Writing result" << std::endl;
  typedef itk::ImageFileWriter<TensorImageType> TensorFileWriterType;
  TensorFileWriterType::Pointer twrit = TensorFileWriterType::New();
  twrit->SetInput(expf->GetOutput());
  twrit->SetFileName(ofile.c_str());
  twrit->Update();
  
}
