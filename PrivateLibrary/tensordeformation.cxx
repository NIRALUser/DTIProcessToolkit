#include "tensordeformation.h"
#include "transforms.h"

#include <itkVectorResampleImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkVectorNearestNeighborInterpolateImageFunction.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

#include "itkDeformationFieldJacobianFilter.h"
#include "itkHFieldToDeformationFieldImageFilter.h"
#include "itkLogEuclideanTensorImageFilter.h"
#include "itkExpEuclideanTensorImageFilter.h"
#include "itkVectorBSplineInterpolateImageFunction.h"
#include "itkTensorRotateImageFilter.h"
#include "itkTensorRotateFromDeformationFieldImageFilter.h"
#include "itkTensorRotateFromDeformationFieldPPDImageFilter.h"

TensorImageType::Pointer createROT(TensorImageType::Pointer timg, 
                                   const std::string &doffile)
{
  RViewTransform<TransformRealType> dof(readDOFFile<TransformRealType>(doffile));
  AffineTransformType::Pointer transform = 
    createITKAffine(dof,
                    timg->GetLargestPossibleRegion().GetSize(),
                    timg->GetSpacing(),
                    timg->GetOrigin());
  vnl_matrix<TransformRealType> R = 
    getInverseRotation(transform);

  typedef itk::TensorRotateImageFilter<TensorImageType,TensorImageType,TransformRealType> TensorRotateType;
  TensorRotateType::Pointer trotfilt = TensorRotateType::New();
  trotfilt->SetRotation(R);
  trotfilt->SetInput(timg);
  trotfilt->Update();

  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType LogTensorImageType;
  LogEuclideanFilter::Pointer logf = LogEuclideanFilter::New();
  logf->SetInput(trotfilt->GetOutput()); 
  logf->Update();

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

  typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
  ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
  expf->SetInput(resampler->GetOutput());
  expf->Update();

  return expf->GetOutput();
    
}

TensorImageType::Pointer createWarp(TensorImageType::Pointer timg,
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

    fsrotate->SetInput1(timg);
    
    fsrotate->SetInput2(jacobian->GetOutput());

    rotate = fsrotate;
    }
  else if(reorientationtype == PreservationPrincipalDirection)
    {
    TensorPPDRotateImageFilterType::Pointer fsrotate;
    fsrotate = TensorPPDRotateImageFilterType::New();

    fsrotate->SetInput1(timg);
    
    fsrotate->SetInput2(jacobian->GetOutput());
    rotate = fsrotate;
    }

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

  typedef itk::ExpEuclideanTensorImageFilter<RealType> ExpEuclideanFilter;
  ExpEuclideanFilter::Pointer expf = ExpEuclideanFilter::New();
  expf->SetInput(warp->GetOutput());
  expf->Update();

  return expf->GetOutput();
  
}

