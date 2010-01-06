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
// ITK includes
#include <itkVersion.h>
#include <itkIndex.h>
#include "fiberio.h"
#include "dtitypes.h"
#include "pomacros.h"
#include "fiberstatsCLP.h"

int main(int argc, char* argv[])
{
#if 0
  namespace po = boost::program_options;
  using namespace boost::lambda;

  // Read program options/configuration
  po::options_description config("Usage: fiberstats input-fiber [options]");
  config.add_options()
    ("help,h", "produce this help message")
    ("verbose,v", "produces verbose output")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("fiber-file", po::value<std::string>(), "DTI fiber file")
    ;

  po::options_description all;
  all.add(config).add(hidden);

  po::positional_options_description p;
  p.add("fiber-file",1);

  po::variables_map vm;

  try
  {
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);     
  } 
  catch (const po::error &e)
  {
    std::cout << "Parse error: " << std::endl;
    std::cout << config << std::endl;
    return EXIT_FAILURE;
  }
#endif
  PARSE_ARGS;
  // End option reading configuration

  // Display help if asked or program improperly called
  if(fiberFile == "")
  {
    std::cout << "No DTI Fiber file given" << std::endl;
      return EXIT_FAILURE;
  }

  const bool VERBOSE = verbose;
  GroupType::Pointer group = readFiberFile(fiberFile);

  verboseMessage("Getting spacing");

  // Get Spacing and offset from group
  const double* spacing = group->GetSpacing();

  const itk::Vector<double, 3> sooffset = 
    group->GetObjectToParentTransform()->GetOffset();

  typedef itk::Index<3> IndexType;
  typedef itk::Functor::IndexLexicographicCompare<3> IndexCompare;
  typedef std::set<IndexType, IndexCompare> VoxelSet;
  typedef std::list<float> MeasureSample;
  typedef std::map<std::string, MeasureSample> SampleMap;

  VoxelSet seenvoxels;
  SampleMap bundlestats;
  bundlestats["fa"] = MeasureSample();
  bundlestats["md"] = MeasureSample();
  bundlestats["fro"] = MeasureSample();

  // For each fiber
  ChildrenListType* children = group->GetChildren(0);
  ChildrenListType::iterator it;
  for(it = children->begin(); it != children->end(); it++)
  {
    DTIPointListType pointlist = 
      dynamic_cast<DTITubeType*>((*it).GetPointer())->GetPoints();
    DTITubeType::Pointer newtube = DTITubeType::New();   

    // For each point along the fiber
    for(DTIPointListType::iterator pit = pointlist.begin();
        pit != pointlist.end(); ++pit)
    {
      typedef DTIPointType::PointType PointType;
      PointType p = pit->GetPosition();

      IndexType i;
      i[0] = static_cast<long int>(round(p[0]));
      i[1] = static_cast<long int>(round(p[1]));
      i[2] = static_cast<long int>(round(p[2]));
      
      seenvoxels.insert(i);

      typedef DTIPointType::FieldListType FieldList;
      const FieldList & fl = pit->GetFields();
      for(FieldList::const_iterator it = fl.begin();
          it != fl.end(); ++it)
      {
        if(bundlestats.count(it->first))
        {
          bundlestats[it->first].push_back(it->second);
        }
      }

    } // end point loop

  } // end fiber loop


  double voxelsize = spacing[0] * spacing[1] * spacing[2];
  std::cout << "Volume (mm^3): " << seenvoxels.size() * voxelsize << std::endl;

  //std::cout << "Measure statistics" << std::endl;
  for(SampleMap::const_iterator it = bundlestats.begin();
      it != bundlestats.end(); ++it)
  {
    const std::string statname = it->first;
    
    double mean = std::accumulate(it->second.begin(), it->second.end(), 0.0) / it->second.size();
    std::cout << statname << " mean: " << mean << std::endl;
    double var = 0.0; //= std::accumulate(it->second.begin(), it->second.end(), 0.0,
      //                                 (_1 - mean)*(_1 - mean)) / (it->second.size() - 1);
    for(MeasureSample::const_iterator it2 = it->second.begin();
        it2 != it->second.end(); it2++)
      {
      double minusMean((*it2) - mean);
      var += (minusMean * minusMean) / (it->second.size() - 1);
      }
    std::cout << statname << " std: " << std::sqrt(var) << std::endl;
  }

  delete children;
  return EXIT_SUCCESS;
}
