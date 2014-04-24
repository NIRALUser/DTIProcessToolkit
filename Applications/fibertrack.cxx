/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.4 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "fiberio.h"
#include "pomacros.h"
#include "itkImageToDTIStreamlineTractographyFilter.h"

#include <itkDiffusionTensor3D.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectWriter.h>
#include <itkMetaDataObject.h>

#include <iostream>
#include <string>
#include <cmath>

#include "fibertrackCLP.h"

enum IntegrationType { Euler, Midpoint, RK4 };
int main(int argc, char* argv[])
{
  typedef itk::DiffusionTensor3D<double> DiffusionTensor;
  typedef itk::Image<DiffusionTensor, 3> TensorImage;
  typedef itk::Image<unsigned short, 3>  LabelImage;
  typedef itk::GroupSpatialObject<3>     FiberBundle;

  typedef itk::ImageFileReader<TensorImage> TensorImageReader;
  typedef itk::ImageFileReader<LabelImage>  LabelImageReader;

  typedef itk::ImageToDTIStreamlineTractographyFilter<TensorImage, LabelImage, FiberBundle> TractographyFilter;

  PARSE_ARGS;

  if( inputTensor == "" || inputROI == "" || outputFiberFile == "" )
    {
    std::cerr << "Tensor image and roi image needs to be specified." << std::endl;
    return EXIT_FAILURE;
    }
  TensorImageReader::Pointer tensorreader = TensorImageReader::New();
  LabelImageReader::Pointer  labelreader  = LabelImageReader::New();

  tensorreader->SetFileName(inputTensor);
  labelreader->SetFileName(inputROI);

  try
    {
    tensorreader->Update();
    labelreader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  if( verbose )
    {
    tensorreader->GetOutput()->Print(std::cout);
    }

  // Sanity check the ROI and tensor image as they must be consistent
  // for the filter to work correctly
  requireequal( (tensorreader->GetOutput()->GetSpacing() == labelreader->GetOutput()->GetSpacing() ),
                "Image Spacings", force);
  requireequal( (tensorreader->GetOutput()->GetLargestPossibleRegion() ==
                 labelreader->GetOutput()->GetLargestPossibleRegion() ),
                "Image Sizes", force);
  requireequal( (tensorreader->GetOutput()->GetOrigin() == labelreader->GetOutput()->GetOrigin() ),
                "Image Origins", force);
  requireequal( (tensorreader->GetOutput()->GetDirection() == labelreader->GetOutput()->GetDirection() ),
                "Image Orientations", force);

  TractographyFilter::Pointer fibertracker = TractographyFilter::New();
  if( reallyVerbose )
    {
    fibertracker->DebugOn();
    }
  if( wholeBrain )
    {
    fibertracker->WholeBrainOn();
    }
  fibertracker->SetTensorImage(tensorreader->GetOutput() );
  fibertracker->SetROIImage(labelreader->GetOutput() );
  fibertracker->SetSourceLabel(sourceLabel);
  fibertracker->SetTargetLabel(targetLabel);
  fibertracker->SetForbiddenLabel(forbiddenLabel);
  fibertracker->SetMaximumAngleChange(maxAngle);
  fibertracker->SetMinimumFractionalAnisotropy(minFa);
  fibertracker->SetStepSize(stepSize);
  fibertracker->Update();

  try
    {
    writeFiberFile(outputFiberFile, fibertracker->GetOutput() );
    }
  catch( itk::ExceptionObject e )
    {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
