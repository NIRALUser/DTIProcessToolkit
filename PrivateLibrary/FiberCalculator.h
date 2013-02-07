#ifndef FIBERCALCULATOR_H
#define FIBERCALCULATOR_H

#include <list>
#include <boost/shared_ptr.hpp>

#include <itkDTITubeSpatialObjectPoint.h>

#include "dtitypes.h"

// Base class for operations in FiberCalculator
// The purpose of this class is to compute a new DTITubePoint from an
// old point and other relevant data.
// The operator () should be overloaded by subclasses to implement the
// modification and results should be returned by value.
class DTIPointModifier
{
public:
  virtual ~DTIPointModifier() {}

  virtual itk::DTITubeSpatialObjectPoint<3> ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<3>& oldpoint) const = 0;
};

// Point modifier to update the position of a point given a warp field
class DTIPointWarper: public DTIPointModifier
{
private:
  typedef itk::VectorInterpolateImageFunction<DeformationImageType> WarpInterpolateType;
public:
  explicit DTIPointWarper(WarpInterpolateType::Pointer wimage) : m_WarpInterpolate(wimage),
                                                                 m_Spacing(wimage->GetInputImage()->GetSpacing()) {}
  virtual ~DTIPointWarper() {}

  virtual itk::DTITubeSpatialObjectPoint<3> ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<3>& oldpoint) const;

private:
  WarpInterpolateType::Pointer m_WarpInterpolate;
  DeformationImageType::SpacingType m_Spacing;
};

// // Point modifier to update the statistics of a point from a tensor
// // field.  The current position of the point is used.
// class DTIPointStatisticsGatherer: public DTIPointModifier
// {
// public:
//   explicit DTIPointStatisticsGatherer(TensorInterpolateType::Pointer timage);
//   virtual ~DTIPointStatisticsGatherer() {}

//   virtual itk::DTITubeSpatialObjectPoint<3> ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<3>& oldpoint) const;
// };

// Point modifier to clear the attribute data from a point and reset
// the tensor to identity.  This is intended to be used before
// gathering new attribute data from an image.
class DTIPointClearData : public DTIPointModifier
{
public:
  virtual ~DTIPointClearData() {}

  virtual itk::DTITubeSpatialObjectPoint<3> ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<3>& oldpoint) const;
};

// Class to perform modification to a fiber bundle based on a set of
// operations.
// This class is designed to implement a chain of responsibility where
// points are based through the operation chain.
class FiberCalculator
{
public:
  FiberCalculator(): m_OldGroup(NULL), m_NewGroup(NULL), m_PointOperationChain() {}
  explicit FiberCalculator(GroupType::Pointer basegroup): m_OldGroup(basegroup), m_NewGroup(basegroup), m_PointOperationChain() {}

  void AddOperation(boost::shared_ptr<DTIPointModifier> operation);
  void ClearOperations();
  void Revert();
  void Update();

private:
  GroupType::Pointer m_OldGroup;
  GroupType::Pointer m_NewGroup;

  std::list<boost::shared_ptr<DTIPointModifier> > m_PointOperationChain;
};

#endif
