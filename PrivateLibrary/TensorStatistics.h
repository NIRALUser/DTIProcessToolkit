// -*- Mode: C++ -*-
/*=============================================================================
  File: TensorStatistics.h
  Author: Tom Fletcher

  TensorStatistics computes statistics of diffusion tensors, including mean,
  covariance, and principal geodesic analysis (PGA). (PGA is the generalization
  of PCA to manifold geometries). The class takes a pointer to a TensorGeometry
  class, which defines how the statistics are computed. In other words, the
  statistics may be linear or nonlinear based on the definition of the
  geometry. For details on the statistics computations, see the reference

  P.T. Fletcher, S. Joshi. Principal Geodesic Analysis on Symmetric Spaces:
  Statistics of Diffusion Tensors. CVAMIA Workshop part of ECCV 2004,
  LNCS 3117, Springer-Verlag, pp. 87-98, 2004.
  http://www.cs.unc.edu/~fletcher/papers/FletcherCVAMIA04_DTStats.pdf

=============================================================================*/

#ifndef __TensorStatistics_h
#define __TensorStatistics_h

#include <itkVectorContainer.h>
#include "TensorGeometry.h"
#include <itkDiffusionTensor3D.h>

template<class T, unsigned int dimension = 3>
class TensorStatistics
{
public:
  typedef itk::DiffusionTensor3D<T> TensorType;
  typedef itk::DiffusionTensor3D<T> TangentType;
  typedef itk::VectorContainer<unsigned int, TensorType> TensorListType;
  typedef typename TensorListType::Pointer TensorListPointerType;

  typedef itk::VectorContainer<unsigned int, T> ScalarListType;
  typedef typename ScalarListType::Pointer ScalarListPointerType;

  typedef itk::SymmetricSecondRankTensor<T, dimension * (dimension + 1) / 2>
    CovarianceType;
  typedef typename CovarianceType::EigenValuesArrayType PGAVariancesArrayType;
  typedef typename CovarianceType::EigenVectorsMatrixType PGAVectorsMatrixType;

  TensorStatistics(TensorGeometry<T, dimension> * _tensGeometry,
                   const T & _stepSize = 1.0)
  {
    tensGeometry = _tensGeometry;
    stepSize = _stepSize;
  }

  void ComputeMean(const TensorListPointerType tensorList,
                   TensorType & mean) const;

  void ComputeWeightedAve(const ScalarListPointerType weightList,
                          const TensorListPointerType tensorList,
                          TensorType & weightedAve) const;

  void ComputeMeanAndCovariance(const TensorListPointerType tensorList,
                                TensorType & mean,
                                CovarianceType & covariance) const;

  void ComputeMeanAndPGA(const TensorListPointerType tensorList,
                         TensorType & mean,
                         PGAVariancesArrayType & pgaVariances,
                         PGAVectorsMatrixType & pgaVectorsMatrix) const;

  // Generates random tensor from isotropic Gaussian density.
  // Depends on the tensor geometry being used.
  TensorType RandomGaussianTensor(const TensorType & mean, T variance) const;

private:
  TensorGeometry<T, dimension> * tensGeometry;
  T stepSize;

  static const T EPSILON = 1.0e-12;
};

#include "TensorStatistics.txx"

#endif
