#include "transforms.h"

#include <itkAffineTransform.h>
#include <itkTransformFileReader.h>
#include <vnl/algo/vnl_svd.h>

#include <iostream>
#include <fstream>

template<class Precision>
RViewTransform<Precision> readDOFFile(const std::string &doffile)
{
  RViewTransform<Precision> dof;
  
  std::ifstream dofstream(doffile.c_str());
  float nil;
  std::string junk;
  dofstream >> junk >> dof.ndofs;
  
  for (unsigned int i = 0; i < 3; ++i)
    dofstream >> nil >> nil >> dof.tx[i];

  for (unsigned int i = 0; i < 3; ++i)
    dofstream >> nil >> nil >> dof.ea[i];

  if(dof.ndofs > 6)
    {
    for (unsigned int i = 0; i < 3; ++i)
      {
      dofstream >> nil >> nil >> dof.sc[i];
      }
    }
  if(dof.ndofs > 9)
    {
    for (unsigned int i = 0; i < 6; ++i)
      dofstream >> nil >> nil >> dof.skew[i];
    }

  return dof;
}

//Read the new format of dof files, it has to be converted first into ASCII format using dof2mat
template<class Precision>
newRViewTransform<Precision> readDOF2MATFile(const std::string &dof2matfile)
{                            
  newRViewTransform<Precision> dofmatrix;
  
  std::ifstream dofstream(dof2matfile.c_str());
  std::string junk;
  //Skip the first line which has 4 characters
  for(unsigned int i = 0 ; i < 4 ; i++)  
    dofstream >> junk;
  
  for (unsigned int i = 0; i < 4; ++i)
    dofstream >> dofmatrix.transfomat(i,0) >> dofmatrix.transfomat(i,1) >> dofmatrix.transfomat(i,2) >> dofmatrix.transfomat(i,3);

  return dofmatrix;
}

template<class Precision, unsigned int ImageDimension>
typename itk::AffineTransform<Precision,
                              ImageDimension>::Pointer
readITKAffine(const std::string &doffile)
{
  typedef itk::AffineTransform<Precision, ImageDimension> ReturnTransformType;

  typedef itk::AffineTransform<float, ImageDimension> FloatTransformType;
  typedef itk::AffineTransform<double, ImageDimension> DoubleTransformType;

  typedef itk::TransformFileReader TransformReader;
  TransformReader::Pointer treader = TransformReader::New();
  treader->SetFileName(doffile);
  treader->Update();
   
  typename ReturnTransformType::Pointer rttransform = ReturnTransformType::New();

  // If we have a float transform
  typename FloatTransformType::Pointer flttransform = 
    dynamic_cast<FloatTransformType*>( treader->GetTransformList()->front().GetPointer() );
  if(!flttransform.IsNull())
  {
    rttransform->SetParameters(flttransform->GetParameters());
    rttransform->SetFixedParameters(flttransform->GetFixedParameters());

    return rttransform;
  }

  // If we have a double transform
  typename DoubleTransformType::Pointer dbltransform = 
    dynamic_cast<DoubleTransformType*>( treader->GetTransformList()->front().GetPointer() );
  if(!dbltransform.IsNull())
  {
    rttransform->SetParameters(dbltransform->GetParameters());
    rttransform->SetFixedParameters(dbltransform->GetFixedParameters());

    return rttransform;
  }
  
  throw itk::ExceptionObject("Invalid transform type");
}


template<class DOFType, class ImageSizeType, class ImageSpacingType, class ImagePointType>
typename itk::AffineTransform<typename DOFType::Precision,
                              ImageSpacingType::Dimension>::Pointer 
createITKAffine(const DOFType & dof,
                const ImageSizeType& size,
                const ImageSpacingType& spacing,
                const ImagePointType& origin)
{
  typedef typename DOFType::Precision Precision;
  typedef typename itk::AffineTransform<Precision,ImageSpacingType::Dimension> AffineTransformType;

  double sinx = sin( dof.ea[0] * M_PI / 180.0);
  double cosx = cos( dof.ea[0] * M_PI / 180.0);
  double siny = sin( dof.ea[1] * M_PI / 180.0);
  double cosy = cos( dof.ea[1] * M_PI / 180.0);
  double sinz = sin( dof.ea[2] * M_PI / 180.0);
  double cosz = cos( dof.ea[2] * M_PI / 180.0);

  vnl_matrix<Precision> rot(3,3);
  rot(0,0) = cosy*cosz;
  rot(0,1) = cosy*sinz;
  rot(0,2) = -siny;
  rot(1,0) = (sinx*siny*cosz-cosx*sinz);
  rot(1,1) = (sinx*siny*sinz+cosx*cosz);
  rot(1,2) = sinx*cosy;
  rot(2,0) = (cosx*siny*cosz+sinx*sinz);
  rot(2,1) = (cosx*siny*sinz-sinx*cosz);
  rot(2,2) = cosx*cosy;

  vnl_matrix<Precision> affine(4,4,0);
  affine(0,0) = cosy*cosz;
  affine(0,1) = cosy*sinz;
  affine(0,2) = -siny;
  affine(0,3) = dof.tx[0];
  affine(1,0) = (sinx*siny*cosz-cosx*sinz);
  affine(1,1) = (sinx*siny*sinz+cosx*cosz);
  affine(1,2) = sinx*cosy;
  affine(1,3) = dof.tx[1];
  affine(2,0) = (cosx*siny*cosz+sinx*sinz);
  affine(2,1) = (cosx*siny*sinz-sinx*cosz);
  affine(2,2) = cosx*cosy;
  affine(2,3) = dof.tx[2];
  affine(3,3) = 1.0;

  if(dof.ndofs > 6)
    {
    affine(0,0) = affine(0,0)*dof.sc[0];
    affine(0,1) = affine(0,1)*dof.sc[1];
    affine(0,2) = affine(0,2)*dof.sc[2];
    affine(1,0) = affine(1,0)*dof.sc[0];
    affine(1,1) = affine(1,1)*dof.sc[1];
    affine(1,2) = affine(1,2)*dof.sc[2];
    affine(2,0) = affine(2,0)*dof.sc[0];
    affine(2,1) = affine(2,1)*dof.sc[1];
    affine(2,2) = affine(2,2)*dof.sc[2];
    
    }

  if(dof.ndofs > 9)
    {
    vnl_matrix<Precision> skewx(4,4,0);
    skewx.set_identity();
    skewx(2, 1)= tan(dof.skew[3]*(M_PI/180.0));
    skewx(1, 2)= tan(dof.skew[2]*(M_PI/180.0)); 
    
    vnl_matrix<Precision> skewy(4,4,0);
    skewy.set_identity();
    skewy(2, 0)= tan(dof.skew[4]*(M_PI/180.0));
    skewy(0, 2)= tan(dof.skew[5]*(M_PI/180.0));
    
    vnl_matrix<Precision> skewz(4,4,0);
    skewz.set_identity();
    skewz(1, 0)= tan(dof.skew[0]*(M_PI/180.0));  
    skewz(0, 1)= tan(dof.skew[1]*(M_PI/180.0)); 

    affine *= skewx * skewy * skewz;

    }

  // itk needs the inverse transform
  vnl_matrix<Precision> affine3(vnl_matrix_inverse<Precision>(affine.extract(3,3)));

  // Setup the itk affine transform 
  typename AffineTransformType::Pointer itktransform = AffineTransformType::New();
  typedef typename AffineTransformType::MatrixType MatrixType;
  typedef typename AffineTransformType::TranslationType TranslationType;
  typedef typename AffineTransformType::CenterType CenterType;
  typedef typename AffineTransformType::OffsetType OffsetType;
  
  MatrixType aff3itk;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      aff3itk(i,j) = affine3(i,j);

  itktransform->SetMatrix(aff3itk);

  // (Extent - 1)/2 * spacing
  CenterType itkcenter;
  itkcenter[0] = (size[0]-1)/2.0 * spacing[0] + origin[0];
  itkcenter[1] = (size[1]-1)/2.0 * spacing[1] + origin[1];
  itkcenter[2] = (size[2]-1)/2.0 * spacing[2] + origin[2];
  itktransform->SetCenter(itkcenter);

  vnl_vector<Precision> tx(3);
  tx[0] = dof.tx[0];
  tx[1] = dof.tx[1];
  tx[2] = dof.tx[2];

  vnl_vector<Precision> rt(3);
  rt = affine3 * tx;

  // the translation needs to be inverted as well
  TranslationType itktranslation;
  itktranslation[0] = -rt[0];
  itktranslation[1] = -rt[1];
  itktranslation[2] = -rt[2];
  itktransform->SetTranslation(itktranslation);  

  return itktransform;
}

template<class DOFType, class ImageSizeType, class ImageSpacingType, class ImagePointType>
typename itk::AffineTransform<typename DOFType::Precision,
                              ImageSpacingType::Dimension>::Pointer 
createnewITKAffine(const DOFType & dof,
                const ImageSizeType& size,
                const ImageSpacingType& spacing,
                const ImagePointType& origin)
{
  typedef typename DOFType::Precision Precision;
  typedef typename itk::AffineTransform<Precision,ImageSpacingType::Dimension> AffineTransformType;

  vnl_matrix<Precision> affine(4,4,0);

  affine(0,0) = dof.transfomat(0,0);
  affine(0,1) = dof.transfomat(0,1);
  affine(0,2) = dof.transfomat(0,2);
  affine(0,3) = dof.transfomat(0,3);
  affine(1,0) = dof.transfomat(1,0);
  affine(1,1) = dof.transfomat(1,1);
  affine(1,2) = dof.transfomat(1,2);
  affine(1,3) = dof.transfomat(1,3);
  affine(2,0) = dof.transfomat(2,0);
  affine(2,1) = dof.transfomat(2,1);
  affine(2,2) = dof.transfomat(2,2);
  affine(2,3) = dof.transfomat(2,3);
  affine(3,3) = dof.transfomat(3,3);


  // itk needs the inverse transform
//  vnl_matrix<Precision> affine3(vnl_matrix_inverse<Precision>(affine.extract(3,3)));

  // Setup the itk affine transform 
  typename AffineTransformType::Pointer itktransform = AffineTransformType::New();
  typedef typename AffineTransformType::MatrixType MatrixType;
  typedef typename AffineTransformType::TranslationType TranslationType;
  typedef typename AffineTransformType::CenterType CenterType;
  typedef typename AffineTransformType::OffsetType OffsetType;
  
  MatrixType aff3itk;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      aff3itk(i,j) = affine(i,j);
  itktransform->SetMatrix(aff3itk);

  // (Extent - 1)/2 * spacing
  CenterType itkcenter;
  itkcenter[0] = (size[0]-1)/2.0 * spacing[0] + origin[0];
  itkcenter[1] = (size[1]-1)/2.0 * spacing[1] + origin[1];
  itkcenter[2] = (size[2]-1)/2.0 * spacing[2] + origin[2];
  itktransform->SetCenter(itkcenter);

  vnl_vector<Precision> tx(3);
  tx[0] = dof.transfomat(0,3);
  tx[1] = dof.transfomat(1,3);
  tx[2] = dof.transfomat(2,3);

  // the translation needs to be inverted as well
  TranslationType itktranslation;
  itktranslation[0] = tx[0];
  itktranslation[1] = tx[1];
  itktranslation[2] = tx[2];
  itktransform->SetTranslation(itktranslation);  

  return itktransform;

}


template <class T>
vnl_matrix<typename T::ObjectType::ScalarType>
getInverseRotation(const T& transform)
{
  typedef typename T::ObjectType::ScalarType Precision;

  // I think we need the inverse for rotating
  vnl_matrix<Precision> inv = vnl_matrix_inverse<Precision>(transform->GetMatrix().GetVnlMatrix());
  vnl_svd<Precision> svd(inv);

  vnl_matrix<Precision> R(svd.U() * svd.V().transpose());

  // R is the rotation including the "rotation" of the skews

  return R;
}
