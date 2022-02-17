// -*- Mode: C++ -*-
/*=============================================================================
  File: MetricSpace.h
  Author: Tom Fletcher

  MetricSpace is an abstract class for defining distances between objects. It's
  only member function is Distance().

=============================================================================*/

#ifndef __MetricSpace_h
#define __MetricSpace_h

#include <cmath>

template <class PointType, class DistanceType = double>
class MetricSpace
{
public:
  virtual ~MetricSpace()
  {
  }

  virtual DistanceType Distance(const PointType & p1, const PointType & p2) = 0;

};

// L2MetricSpace defines the L^2 metric for vector/scalar types.
template <class PointType, class DistanceType = double>
class L2MetricSpace : public MetricSpace<PointType, DistanceType>
{
public:
  virtual ~L2MetricSpace()
  {
  }

  virtual DistanceType Distance(const PointType & p1, const PointType & p2);

};

template <class PointType, class DistanceType>
DistanceType
L2MetricSpace<PointType, DistanceType>::Distance(const PointType & p1,
                                                 const PointType & p2)
{
  return sqrt( (p1 - p2) * (p1 - p2) );
}

#endif
