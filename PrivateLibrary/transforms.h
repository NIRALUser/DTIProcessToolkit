#ifndef _transforms_h__
#define _transforms_h__

#include <itkAffineTransform.h>
#include <itkMatrix.h>
#include <vnl/vnl_matrix.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI


template<class TPrecision>
class RViewTransform
{
public:
  typedef TPrecision Precision;
  Precision tx[3];  // Translation
  Precision ea[3];  // Euler angles
  Precision sc[3];  // Scaling 
  Precision skew[6];// Skews
  unsigned int ndofs;
};

template<class TPrecision>
class newRViewTransform
{
public:
  typedef TPrecision Precision;
  itk::Matrix<Precision,4,4> transfomat;
};

// Reads the parameter's from RView's dof file
template<class Precision>
RViewTransform<Precision>
readDOFFile(const std::string &doffile);

// Reads the parameter's from the new version of RView's dof file
template<class Precision>
newRViewTransform<Precision>
readDOF2MATFile(const std::string &doffile);

template<class Precision, unsigned int ImageDimension>
typename itk::AffineTransform<Precision,
                              ImageDimension>::Pointer
readITKAffine(const std::string &doffile);

// Creates an ITK affine transform from the RView's dof file
template <class DOFType, 
          class ImageSizeType, 
          class ImageSpacingType, 
          class ImagePointType>
typename itk::AffineTransform<typename DOFType::Precision, 
                              ImageSizeType::Dimension>::Pointer 
createITKAffine(const DOFType& dof,
                const ImageSizeType& size,
                const ImageSpacingType& spacing,
                const ImagePointType& origin);

// Creates an ITK affine transform from the NEW RView's dof file
template <class DOFType, 
          class ImageSizeType, 
          class ImageSpacingType, 
          class ImagePointType>
typename itk::AffineTransform<typename DOFType::Precision, 
                              ImageSizeType::Dimension>::Pointer 
createnewITKAffine(const DOFType& dof,
                const ImageSizeType& size,
                const ImageSpacingType& spacing,
                const ImagePointType& origin);

// Returns the inverse of the rotational component
// of an affine transformation
template <class T>
vnl_matrix<typename T::ObjectType::ScalarType>
getInverseRotation(const T &transform);

#include "transforms.txx"

#endif
