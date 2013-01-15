/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2008-08-06 12:51:30 -0400 (Wed, 06 Aug 2008) $
  Version:   $Revision: 597 $
  Author:    Marc Niethammer (mn@cs.unc.edu)
             Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// This program creates a DWI atlas from a set of DWIs and associated
// deformation maps to a common atlas space

// Authors: Sami Benzaid, Casey Goodlett, Marc Niethammer, 2010

// STL includes
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>

// boost includes
#include <boost/program_options/option.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/cmdline.hpp>

#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

// ITK includes
// datastructures
#include <itkImage.h>
#include <itkVector.h>
#include <itkMetaDataObject.h>

// IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Filters
#include <itkTensorFractionalAnisotropyImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

#include <itkNthElementImageAdaptor.h>
#include <itkOtsuThresholdImageCalculator.h>

#include "itkVectorMaskNegatedImageFilter.h"
#include "itkVectorMaskImageFilter.h"
#include "itkDiffusionTensor3DReconstructionLinearImageFilter.h"

#include "itkHFieldToDeformationFieldImageFilter.h"
#include "itkDeformationFieldJacobianFilter.h"

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_inverse.h>

#include "dtitypes.h"
#include "deformationfieldio.h"

#include "dwiAtlasCLP.h"

#include "auxFunctions.h"
#include "sphericalHarmonicsFunctions.h"

#include <itkVariableLengthVector.h>
#include <ctime>

#include "itkDWIAtlasBuilder.h"

const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

namespace po = boost::program_options;

// Define necessary types for images
//typedef double MyRealType;
typedef float MyRealType;

//const unsigned int DIM = 3;

typedef unsigned short DWIPixelType;
typedef itk::VectorImage<DWIPixelType, DIM> VectorImageType;
typedef itk::Image<DWIPixelType, DIM> ScalarImageType;
typedef itk::Image<unsigned char, DIM> MaskImageType;
typedef itk::Image<DWIPixelType, DIM> IntImageType;
typedef itk::DiffusionTensor3DReconstructionLinearImageFilter<DWIPixelType, MyRealType>
DiffusionEstimationFilterType;
typedef DiffusionEstimationFilterType::GradientDirectionContainerType GradientDirectionContainerType;
typedef vnl_vector_fixed<MyRealType, 3> MyGradientType;

typedef itk::DWIAtlasBuilder< MyRealType, DWIPixelType > DWIAtlasBuilderType;

int outputAllSettings(int argc, char* argv[] )
{
  PARSE_ARGS;
  
  std::cout << "Parameter settings:" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "output volume        = " << dwiOutputVolume << std::endl;
  std::cout << "gradient vector file = " << gradientVectorFile << std::endl;
  std::cout << "case file            = " << sCaseFile << std::endl;
  std::cout << "outlierImageName     = " << outlierImageName << std::endl;
  std::cout << "maskImageFileName    = " << maskImageFileName << std::endl;
  std::cout << "h-field              = " << bHField << std::endl;
  std::cout << "bJustDoResampling    = " << bJustDoResampling << std::endl;
  std::cout << "verbose              = " << bVERBOSE << std::endl;
  std::cout << "SHOrder              = " << SHOrder << std::endl;
  std::cout << "lambda               = " << dLambda << std::endl;
  std::cout << "nrOfThreads          = " << nrOfThreads << std::endl;
  std::cout << "interpolationType    = " << interpolationType << std::endl;
  std::cout << "averagingType        = " << averagingType << std::endl;
  std::cout << "bNoLogFit            = " << bNoLogFit << std::endl;
  std::cout << "dScalingFactor       = " << dScalingFactor << std::endl;
  std::cout << "dHuberC              = " << dHuberC << std::endl;
  std::cout << "bDoWeightedLSEst     = " << bDoWeightedLSEst << std::endl;
  std::cout << "dRiceSigma           = " << dRiceSigma << std::endl;
  std::cout << "iLogMinArgumentValue = " << iLogMinArgumentValue << std::endl;
  std::cout << "iWLSIter             = " << iWLSIter << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << std::endl << std::endl << std::endl;

  return 0;

}

int main(int argc, char* argv[])
{

//  unsigned int scale;

  PARSE_ARGS;

  if ( bVERBOSE )
    outputAllSettings( argc, argv );

  DWIAtlasBuilderType::Pointer pAtlasBuilder = DWIAtlasBuilderType::New();

  // Set all the command-line arguments we specified
  
  // the case file is the input; the filter takes care of loading and
  // processing everything
  pAtlasBuilder->SetInput( &sCaseFile );  

  pAtlasBuilder->SetVerbose( bVERBOSE );
  pAtlasBuilder->SetIsHField( bHField );
  pAtlasBuilder->SetJustDoResampling( bJustDoResampling );
  pAtlasBuilder->SetGradientVectorFile( gradientVectorFile );
  pAtlasBuilder->SetSHOrder( SHOrder );
  pAtlasBuilder->SetLambda( dLambda );
  pAtlasBuilder->SetDoWeightedLS( bDoWeightedLSEst );
  pAtlasBuilder->SetRiceSigma( (MyRealType)dRiceSigma );
  pAtlasBuilder->SetHuberC( (MyRealType)dHuberC );
  pAtlasBuilder->SetInterpolationType( interpolationType );
  pAtlasBuilder->SetAveragingType( averagingType );
  pAtlasBuilder->SetNrOfThreads( nrOfThreads );
  pAtlasBuilder->SetNoLogFit( bNoLogFit );
  pAtlasBuilder->SetScalingFactor( (MyRealType)dScalingFactor );
  pAtlasBuilder->SetLogMinArgumentValue( (DWIPixelType)iLogMinArgumentValue );
  pAtlasBuilder->SetNrOfWLSIterations( (unsigned int)iWLSIter );
  pAtlasBuilder->SetMaskImageFileName( maskImageFileName );

  if ( bNoLogFit && bDoWeightedLSEst )
    {
    std::cout << "ERROR: WLS in the original signal domain is not supported. Remove --noLogFit specifier and try again. ABORT." << std::endl;
    return EXIT_FAILURE;
    }

  if ( iLogMinArgumentValue != 1 )
    {
    std::cout << "WARNING: iLogMinArgumentValue should always be one, but it is " << iLogMinArgumentValue << std::endl;
    std::cout << "Unless you really know what you are doing, do NOT proceed." << std::endl;
    }

  if ( bVERBOSE )
    {
    using namespace boost::posix_time;
    ptime now = second_clock::local_time();
    // write out time when the computation started
    std::cout << "Processing started at: " << now << std::endl << std::endl;
    }

  pAtlasBuilder->Update();

  if (bVERBOSE)
    std::cerr << "Writing the final result to: " << dwiOutputVolume << std::endl;

  // Write reconstructed DWI file 
  try
  {
    typedef itk::ImageFileWriter<VectorImageType> VectorFileWriterType;

    VectorFileWriterType::Pointer dwiWriter = VectorFileWriterType::New();
    dwiWriter->SetFileName( dwiOutputVolume.c_str() );
  
    dwiWriter->SetInput( pAtlasBuilder->GetOutput() );
    dwiWriter->SetUseCompression(true);
    dwiWriter->Update();
       
  } 
  catch (itk::ExceptionObject e) 
  {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }  

  if ( outlierImageName.compare("None")!=0 )
    {
    if ( bVERBOSE )
      std::cout << "Writing outlier image to " << outlierImageName << " ... ";

    typedef itk::ImageFileWriter<ScalarImageType> ScalarFileWriterType;

    ScalarFileWriterType::Pointer scalarWriter = ScalarFileWriterType::New();
    scalarWriter->SetFileName( outlierImageName.c_str() );

    scalarWriter->SetInput( pAtlasBuilder->GetOutlierImage() );
    scalarWriter->SetUseCompression(true);
    scalarWriter->Update();

    if ( bVERBOSE )
      std::cout << "done." << std::endl;
    }

  return EXIT_SUCCESS;
}

