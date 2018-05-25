/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.3 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// STL includes
#include <string>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

// ITK includes
#include <itkVersion.h>
#include <itkIndex.h>
#include "fiberio.h"
#include "dtitypes.h"
#include "pomacros.h"
#include "fiberstatsCLP.h"

#if ITK_VERSION_MAJOR >= 5
#include <itkLexicographicCompare.h>
#endif

int main(int argc, char* argv[])
{
  PARSE_ARGS;
  // End option reading configuration
  const bool         VERBOSE = verbose;
  GroupType::Pointer group = readFiberFile(fiberFile);

  verboseMessage("Getting spacing");

  // Get Spacing and offset from group
  const double* spacing = group->GetSpacing();

  typedef itk::Index<3>                              IndexType;
#if ITK_VERSION_MAJOR >= 5
  typedef itk::Functor::LexicographicCompare<itk::Index<3>, itk::Index<3> > IndexCompare;
#else
  typedef itk::Functor::IndexLexicographicCompare<3> IndexCompare;
#endif
  typedef std::set<IndexType, IndexCompare>          VoxelSet;
  typedef std::list<float>                           MeasureSample;
  typedef std::map<std::string, MeasureSample>       SampleMap;
  typedef std::map<IndexType, int, IndexCompare>         VoxelMap;
  typedef VoxelMap::value_type    VoxelMapValue;

  VoxelSet  seenvoxels;
  VoxelMap  visitvoxels;
  SampleMap bundlestats;
  bundlestats["fa"] = MeasureSample();
  bundlestats["md"] = MeasureSample();
  bundlestats["l1"] = MeasureSample();
  bundlestats["rd"] = MeasureSample();

  // For each fiber
  std::vector<double> FiberLengthsVector;
  ChildrenListType*          children = group->GetChildren(0);
  ChildrenListType::iterator it;

  for( it = children->begin(); it != children->end(); it++ )
    {
    DTIPointListType pointlist =
      dynamic_cast<DTITubeType *>( (*it).GetPointer() )->GetPoints();
    DTITubeType::Pointer newtube = DTITubeType::New();
    // For each point along the fiber
    double FiberLength=0;// Added by Adrien Kaiser 04-03-2013
    typedef DTIPointType::PointType PointType;
    PointType Previousp = pointlist.begin()->GetPosition();// Added by Adrien Kaiser 04-03-2013
    for( DTIPointListType::iterator pit = pointlist.begin();
         pit != pointlist.end(); ++pit )
      {
      PointType p = pit->GetPosition();

      // Added by Adrien Kaiser 04-03-2013: Compute length between 2 points
      if(pit != pointlist.begin()) // no previous for the first one
      {
        double length = sqrt( (Previousp[0]-p[0])*(Previousp[0]-p[0]) + (Previousp[1]-p[1])*(Previousp[1]-p[1]) +(Previousp[2]-p[2])*(Previousp[2]-p[2]) );
        FiberLength = FiberLength + length;
      }
      Previousp = p;
      //

      IndexType i;
      i[0] = static_cast<long int>(vnl_math_rnd_halfinttoeven(p[0]) );
      i[1] = static_cast<long int>(vnl_math_rnd_halfinttoeven(p[1]) );
      i[2] = static_cast<long int>(vnl_math_rnd_halfinttoeven(p[2]) );

      seenvoxels.insert(i);
      int curValue = visitvoxels[i];
      if (curValue == 0)
      { 
	visitvoxels[i] = 1;
      } else {
	visitvoxels[i] = curValue + 1;
      }

      typedef DTIPointType::FieldListType FieldList;
      const FieldList & fl = pit->GetFields();
      for( FieldList::const_iterator flit = fl.begin();
           flit != fl.end(); ++flit )
        {
        if( bundlestats.count(flit->first) )
          {
          bundlestats[flit->first].push_back(flit->second);
          }
        }

      } // end point loop
    FiberLengthsVector.push_back(FiberLength);// Added by Adrien Kaiser 04-03-2013

    }   // end fiber loop

  // Added by Adrien Kaiser 04-03-2013: compute average fiber length and quantiles
  std::cout<< FiberLengthsVector.size() <<" fibers found"<<std::endl;

  if( FiberLengthsVector.empty() )
  {
    std::cout<<"This fiber file is empty. ABORT."<<std::endl;
    return EXIT_FAILURE;
  }

  double AverageFiberLength = 0;
  for(unsigned int AFLVecIter=0; AFLVecIter < FiberLengthsVector.size();AFLVecIter++)
  {
     AverageFiberLength = AverageFiberLength + FiberLengthsVector[AFLVecIter];
  }
  AverageFiberLength = AverageFiberLength / FiberLengthsVector.size();
  std::cout<<"Average Fiber Length: "<<AverageFiberLength<<std::endl;

  sort( FiberLengthsVector.begin(), FiberLengthsVector.end() );
  std::cout<<"Minimum Fiber Length: "<< FiberLengthsVector[0] <<std::endl;
  std::cout<<"Maximum Fiber Length: "<< FiberLengthsVector[FiberLengthsVector.size()-1] <<std::endl;
  std::cout<<"75 percentile Fiber Length: "<< FiberLengthsVector[ (int)(0.75*FiberLengthsVector.size()) ] <<std::endl;
  std::cout<<"90 percentile Fiber Length: "<< FiberLengthsVector[ (int)(0.9*FiberLengthsVector.size()) ] <<std::endl;
  double Average75PercFiberLength = 0;
  for(unsigned int AFLVecIter=(int)(0.75*FiberLengthsVector.size()); AFLVecIter < FiberLengthsVector.size();AFLVecIter++)
  {
     Average75PercFiberLength = Average75PercFiberLength + FiberLengthsVector[AFLVecIter];
  }
  Average75PercFiberLength = Average75PercFiberLength / (FiberLengthsVector.size() - (int)(0.75*FiberLengthsVector.size())) ;
  std::cout<<"Average 75 Percentile Fiber Length: "<<Average75PercFiberLength<<std::endl;
  //

  double voxelsize = spacing[0] * spacing[1] * spacing[2];
  std::cout << "Volume (mm^3): " << seenvoxels.size() * voxelsize  << std::endl;

  float weightedSum = 0.0;
  for ( VoxelMap::const_iterator voxelIter = visitvoxels.begin(); voxelIter != visitvoxels.end(); ++voxelIter)
  {
    weightedSum += voxelIter->second;
  }

  //int weightedSum =  std::accumulate(visitvoxels.begin(), visitvoxels.end(), 0.0);
  std::cout << "Density Volume (mm^3): " << weightedSum * voxelsize / FiberLengthsVector.size()  << std::endl;
  std::cout << "Measure statistics: " << bundlestats.size() << std::endl;
  for( SampleMap::const_iterator smit = bundlestats.begin();
       smit != bundlestats.end(); ++smit )
  {
    const std::string statname = smit->first;

    if (smit->second.size() > 0) 
    {
      double mean = std::accumulate(smit->second.begin(), smit->second.end(), 0.0) / smit->second.size();
      std::cout << statname << " mean: " << mean << std::endl;
      // double var = 0.0; // = std::accumulate(smit->second.begin(), smit->second.end(), 0.0,
      //                                 (_1 - mean)*(_1 - mean)) / (smit->second.size() - 1);
      // for( MeasureSample::const_iterator it2 = smit->second.begin();
      // 	   it2 != smit->second.end(); it2++ )
      // 	{
      // 	double minusMean( (*it2) - mean);
      // 	var += (minusMean * minusMean) / (smit->second.size() - 1);
      // 	}
      // std::cout << statname << " std: " << std::sqrt(var) << std::endl;
    } else {
      std::cout << statname << " has no measures " << std::endl;
    }
      
   }

  delete children;
  return EXIT_SUCCESS;
}
