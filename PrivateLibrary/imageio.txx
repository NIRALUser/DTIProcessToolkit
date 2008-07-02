/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2008-07-02 15:54:54 $
  Version:   $Revision: 1.3 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "imageio.h"
#include <itkImageFileWriter.h>

template<typename TImage>
void writeImage(const std::string & filename, 
                typename itk::SmartPointer<TImage> image)
{
  typedef itk::ImageFileWriter<TImage> ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetUseCompression(true);
  writer->SetInput(image);
  writer->SetFileName(filename.c_str());
  writer->Update();

}
