// -*- Mode: C++ -*-
/*=============================================================================
  File: SymmetricSpaceTensorGeometry.h
  Author: Tom Fletcher

  SymmetricSpaceTensorGeometry defines the symmetric space geometry of the
  space of diffusion tensors. For details on this geometry see the reference

  P.T. Fletcher, S. Joshi. Principal Geodesic Analysis on Symmetric Spaces:
  Statistics of Diffusion Tensors. CVAMIA Workshop part of ECCV 2004,
  LNCS 3117, Springer-Verlag, pp. 87-98, 2004.
  http://www.cs.unc.edu/~fletcher/papers/FletcherCVAMIA04_DTStats.pdf

=============================================================================*/

#ifndef __SymmetricSpaceTensorGeometry_h
#define __SymmetricSpaceTensorGeometry_h

#include <vnl/vnl_det.h>
#include <vnl/vnl_trace.h>
#include "TensorGeometry.h"

template<class T, unsigned int dimension=3>
class SymmetricSpaceTensorGeometry : public TensorGeometry<T, dimension>
{
public:
  typedef TensorGeometry<T, dimension> SuperClass;
  typedef typename SuperClass::TensorType TensorType;
  typedef typename SuperClass::TangentType TangentType;

  typedef itk::Matrix<T, dimension> MatrixType;
  typedef typename TensorType::EigenValuesArrayType EigenValuesArrayType;
  typedef typename TensorType::EigenVectorsMatrixType EigenVectorsMatrixType;

  SymmetricSpaceTensorGeometry() {}

  virtual T InnerProduct(const TensorType & base,
                         const TangentType & v,
                         const TangentType & w);

  virtual TensorType ExpMap(const TensorType & base, const TangentType & v);
  virtual TangentType LogMap(const TensorType & base, const TensorType & p);

  TensorType GroupAction(const TensorType & p, const MatrixType & g);

private:
  // Computes eigensystem where the eigenvector matrix is determinant 1.
  void OrientedEigensystem(const TensorType & p,
                           EigenValuesArrayType & eigenValues,
                           EigenVectorsMatrixType & eigenVectors);

  // Perhaps these routines should be part of SymmetricSecondRankTensor?
  void TensorToMatrix(const TensorType & p, MatrixType & m);
  void MatrixToTensor(const MatrixType & m, TensorType & p);
};

#include "SymmetricSpaceTensorGeometry.txx"

#endif
