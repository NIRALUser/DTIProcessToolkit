// -*- Mode: C++ -*-
/*=============================================================================
  File: TensorGeometry.h
  Author: Tom Fletcher

  TensorGeometry is an abstract class for defining the geometry of the space
  of diffusion tensors. The geometry is defined via a Riemannian metric, and
  this leads to geodesics and distances. Note that the geometry does not need
  to be curved, in other words, a linear geometry is also possible in this
  framework.

=============================================================================*/

#ifndef __TensorGeometry_h
#define __TensorGeometry_h

#include <itkSymmetricSecondRankTensor.h>
#include "MetricSpace.h"

template<class T, unsigned int dimension=3>
class TensorGeometry : public MetricSpace<itk::DiffusionTensor3D<T>, T>
{
public:

  typedef itk::DiffusionTensor3D<T> TensorType;
  typedef itk::DiffusionTensor3D<T> TangentType;

  TensorGeometry() {}

  // Defines the Riemannian metric at a base point for two tangent vectors.
  virtual T InnerProduct(const TensorType & base,
                         const TangentType & v,
                         const TangentType & w) = 0;
  virtual T Norm(const TensorType & base, const TangentType & v);
  virtual T NormSquared(const TensorType & base, const TangentType & v);

  // ExpMap gives the Riemannian exponential map (i.e., geodesic segment
  // starting at base with initial velocity v)
  virtual TensorType ExpMap(const TensorType & base, const TangentType & v) = 0;

  // LogMap is the inverse of the exponential map. It returns the initial
  // veloctiy vector for the geodesic segment between base and p.
  virtual TangentType LogMap(const TensorType & base, const TensorType & p) = 0;

  // Geodesic distance between tensors a and b.
  virtual T Distance(const TensorType & a, const TensorType & b);
};

template<class T, unsigned int dimension>
T TensorGeometry<T, dimension>::Norm(const TensorType & base,
                                     const TangentType & v)
{
  return sqrt(NormSquared(base, v));
}

template<class T, unsigned int dimension>
T TensorGeometry<T, dimension>::NormSquared(const TensorType & base,
                                            const TangentType & v)
{
  return InnerProduct(base, v, v);
}

template<class T, unsigned int dimension>
T TensorGeometry<T, dimension>::Distance(const TensorType & a,
                                         const TensorType & b)
{
  return Norm(a, LogMap(a, b));
}

#endif
