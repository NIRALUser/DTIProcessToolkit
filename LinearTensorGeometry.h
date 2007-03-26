// -*- Mode: C++ -*-
/*=============================================================================
  File: LinearTensorGeometry.h
  Author: Tom Fletcher

  LinearTensorGeometry defines a linear geometry on the space of diffusion
  tensors. The inner product is the standard Frobenius inner product between
  matrices.

=============================================================================*/

#ifndef __LinearTensorGeometry_h
#define __LinearTensorGeometry_h

#include "TensorGeometry.h"

template<class T, unsigned int dimension=3>
class LinearTensorGeometry : public TensorGeometry<T, dimension>
{
public:
  typedef TensorGeometry<T, dimension> SuperClass;
  typedef typename SuperClass::TensorType TensorType;
  typedef typename SuperClass::TangentType TangentType;

  LinearTensorGeometry() {}

  virtual T InnerProduct(const TensorType & base, const TangentType & v,
                         const TangentType & w);

  virtual TensorType ExpMap(const TensorType & base, const TangentType & v);
  virtual TangentType LogMap(const TensorType & base, const TensorType & p);
};

template<class T, unsigned int dimension>
T LinearTensorGeometry<T, dimension>::InnerProduct(const TensorType & base,
                                                   const TangentType & v,
                                                   const TangentType & w)
{
  T prod;
  int i, j, diagIndex;

  // Diagonal elements are weighted by one, off-diagonal by two
  diagIndex = 0;
  j = dimension;
  prod = 0.0;
  for(i = 0; i < v.Size(); i++)
  {
    if(i == diagIndex)
    {
      prod += v[i] * w[i];
      diagIndex += j;
      j--;
    }
    else
      prod += 2.0 * v[i] * w[i];
  }

  return prod;
}

template<class T, unsigned int dimension>
typename LinearTensorGeometry<T, dimension>::TensorType
LinearTensorGeometry<T, dimension>::ExpMap(const TensorType & base,
                                           const TangentType & v)
{
  return (base + v);
}

template<class T, unsigned int dimension>
typename LinearTensorGeometry<T, dimension>::TangentType
LinearTensorGeometry<T, dimension>::LogMap(const TensorType & base,
                                           const TensorType & p)
{
  return (p - base);
}

#endif
