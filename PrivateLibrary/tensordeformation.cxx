#include "tensordeformation.h"
#include "transforms.h"

#include <itkVectorResampleImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkVectorNearestNeighborInterpolateImageFunction.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkTransformFileReader.h>
#include <itkTransformBase.h>

#include "itkDeformationFieldJacobianFilter.h"
#include "itkHFieldToDeformationFieldImageFilter.h"
#include "itkLogEuclideanTensorImageFilter.h"
#include "itkExpEuclideanTensorImageFilter.h"
#include "itkVectorBSplineInterpolateImageFunction.h"
#include "itkTensorRotateImageFilter.h"
#include "itkTensorRotateFromDeformationFieldImageFilter.h"
#include "itkTensorRotateFromDeformationFieldPPDImageFilter.h"

TensorImageType::Pointer createROT(TensorImageType::Pointer timg, 
                                   const std::string &doffile,
				   int doffiletype)
{
  //Depending on which kind of input is given for the transformation we adapt the readers:
  // doffiletype = 0 -> Old dof file, = 1 -> New dof file (dof2mat), = 2 -> itk format.

  AffineTransformType::Pointer transform;
  //Deal with DOF files
  if(doffiletype == 0)
  {
    RViewTransform<TransformRealType> dof(readDOFFile<TransformRealType>(doffile));
    //AffineTransformType::Pointer transform = 
    AffineTransformType::Pointer transform_tmp = 
      createITKAffine(dof,
		      timg->GetLargestPossibleRegion().GetSize(),
		      timg->GetSpacing(),
		      timg->GetOrigin());
    transform = transform_tmp;
  }
  if(doffiletype == 1)
  {
    newRViewTransform<TransformRealType> dof(readDOF2MATFile<TransformRealType>(doffile));
    AffineTransformType::Pointer transform_tmp = 
      createnewITKAffine(dof,
		      timg->GetLargestPossibleRegion().GetSize(),
		      timg->GetSpacing(),
		      timg->GetOrigin());
    transform = transform_tmp;
  }
  //Deal with itk affine files
  else if(doffiletype == 2)
  {
    //Create a temporary transform filter before copying it to transform
    AffineTransformType::Pointer transform_tmp = AffineTransformType::New();

    typedef itk::TransformFileReader TransformationReader;
    typedef itk::TransformBase TransformBaseType;
    TransformationReader::Pointer treader = TransformationReader::New();

    treader->SetFileName(doffile);
    try
    {
      treader->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << e << std::endl;
    }
    
    TransformBaseType::Pointer basetransform = treader->GetTransformList()->front();
   
    //Fill out the Affine matrix
    AffineTransformType::MatrixType aff3itk;
    int x = 0;
    int y = 0;
    for(unsigned int i = 0 ; i < 9 ; i++)
    {
      aff3itk(x,y) = basetransform->GetParameters()(i);
      y++;
      if(y == 3) 
      {
	y = 0;
	x++;
      }
    }
    transform_tmp->SetMatrix(aff3itk);

    //Set the translation values
    AffineTransformType::TranslationType itktranslation;
    itktranslation[0] = basetransform->GetParameters()(9);
    itktranslation[1] = basetransform->GetParameters()(10);
    itktranslation[2] = basetransform->GetParameters()(11);
    transform_tmp->SetTranslation(itktranslation);

    //Get the fixed parameters (center of rotation)
    AffineTransformType::CenterType itkcenter;
    itkcenter[0] = basetransform->GetFixedParameters()(0);
    itkcenter[1] = basetransform->GetFixedParameters()(1);
    itkcenter[2] = basetransform->GetFixedParameters()(2);
    transform_tmp->SetCenter(itkcenter);

    //Copy the temporary transform in the transform filter
    transform = transform_tmp;
  }

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
                                    DeformationImageType::Pointer forward,
                                    TensorReorientationType reorientationtype,
                                    InterpolationType interpolationtype)
{
  // Compute jacobian of inverse deformation field
  typedef itk::DeformationFieldJacobianFilter<DeformationImageType, RealType> JacobianFilterType;
  typedef JacobianFilterType::OutputImageType JacobianImageType;
  JacobianFilterType::Pointer jacobian = JacobianFilterType::New();
  //jacobian->SetInput(inverse);
  jacobian->SetInput(forward);

  typedef itk::LogEuclideanTensorImageFilter<RealType> LogEuclideanFilter;
  typedef LogEuclideanFilter::OutputImageType LogTensorImageType;
  LogEuclideanFilter::Pointer logf = LogEuclideanFilter::New();
  logf->SetInput(timg); 
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
  warp->SetDeformationField(forward);
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

  //return expf->GetOutput();
  
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

    fsrotate->SetInput1(expf->GetOutput());
    
    fsrotate->SetInput2(jacobian->GetOutput());

    rotate = fsrotate;
    }
  else if(reorientationtype == PreservationPrincipalDirection)
    {
    TensorPPDRotateImageFilterType::Pointer ppdrotate;
    ppdrotate = TensorPPDRotateImageFilterType::New();

    ppdrotate->SetInput1(expf->GetOutput());
    
    ppdrotate->SetInput2(jacobian->GetOutput());
    rotate = ppdrotate;
    }
  rotate->Update();
  return rotate->GetOutput();
}

