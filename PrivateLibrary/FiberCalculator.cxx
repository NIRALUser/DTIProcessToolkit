#include "FiberCalculator.h"

itk::DTITubeSpatialObjectPoint<DIM>
DTIPointWarper::ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<DIM> & oldpoint) const
{
  itk::DTITubeSpatialObjectPoint<DIM> newpoint(oldpoint);

  typedef itk::DTITubeSpatialObjectPoint<DIM>::PointType PointType;
  PointType p = oldpoint->GetPosition();
  typedef WarpInterpolateType::ContinuousIndexType ContinuousIndexType;
  ContinuousIndexType ci;

  std::copy(p, p+DIM, ci);

  DeformationPixelType displacement(m_WarpInterpolate->EvaluateAtContinuousIndex(ci).GetDataPointer());

  for(unsigned int i = 0; i < 3; ++i)
    {
    p[i] = p[i] + displacement[i] / m_Spacing[i];
    }
  return newpoint;
}

// itk::DTITubeSpatialObjectPoint<DIM>
// DTIPointStatisticsGatherer::ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<DIM> & oldpoint) const
// {
  
  
// }


itk::DTITubeSpatialObjectPoint<DIM>
DTIPointClearData::ComputeNewPoint(const itk::DTITubeSpatialObjectPoint<DIM> & oldpoint) const
{
  itk::DTITubeSpatialObjectPoint<DIM> newpoint;
  newpoint.SetPosition(oldpoint->GetPosition());
  return newpoint;
}
