#include <itkDiffusionTensor3D.h>
#include <itkVector.h>

template <class TPrecision>
class TensorWithDerivatives
{
public:
  typedef itk::DiffusionTensor3D<TPrecision> TensorType;
  typedef itk::Vector<itk::Vector<TPrecision, 3>, 6> TensorSpatialGradientType;
  typedef itk::Vector<itk::SymmetricSecondRankTensor<TPrecision, 3>, 6>
                                                     TensorSpatialHessianType;
  TensorType D;
  TensorSpatialGradientType gradD;
  TensorSpatialGradientType nablaD;

}
