#ifndef _transforms_h__
#define _transforms_h__

#include <itkAffineTransform.h>
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

// template<class Precision>
// boost::shared_ptr<RViewTransform<Precision> > readDOFFile(const std::string &doffile);

// template <class ImageSizeType, class ImageSpacingType, class Precision>
// typename itk::AffineTransform<Precision,ImageSizeType::Dimension>::Pointer 
// createITKAffine(const boost::shared_ptr<RViewTransform>& dof,
//                 const ImageSizeType& size,
//                 const ImageSpacingType& spacing);

template <class T>
vnl_matrix<typename T::ObjectType::ScalarType>
getInverseRotation(const T &transform);

#include "transforms.txx"

#endif
