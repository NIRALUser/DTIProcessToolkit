/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>

#include <itkImageFileWriter.h>

#include "deformationfieldoperations.h"
#include "dtitypes.h"

int main(int argc, char* argv[])
{
  if( argc != 3)
    {
    std::cerr << "Usage: " << argv[0] << " <hfield> <displacementfield>" << std::endl;
    return EXIT_FAILURE;
    }
  std::string infile(argv[1]);
  std::string outfile(argv[2]);

  DeformationImageType::Pointer deffield = readDeformationField(infile, HField);

  typedef itk::ImageFileWriter<DeformationImageType> DeformationFileWriter;
  DeformationFileWriter::Pointer defwriter = DeformationFileWriter::New();
  defwriter->SetInput(deffield);
  defwriter->SetFileName(outfile);
  defwriter->Update();

  return EXIT_SUCCESS;
}
