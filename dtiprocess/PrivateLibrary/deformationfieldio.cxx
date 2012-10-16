#include "deformationfieldio.h"
#include <string>
#include <itkImageFileReader.h>
#include "itkHFieldToDeformationFieldImageFilter.h"

DeformationImageType::Pointer readDeformationField(std::string warpfile, DeformationFieldType dft)
{
  typedef itk::ImageFileReader<DeformationImageType> DeformationImageReader;

  DeformationImageReader::Pointer defreader = DeformationImageReader::New();
  defreader->SetFileName(warpfile.c_str());

  if(dft == HField)
    {

    typedef itk::HFieldToDeformationFieldImageFilter<DeformationImageType> DeformationConvertType;
    DeformationConvertType::Pointer defconv = DeformationConvertType::New();
    defconv->SetInput(defreader->GetOutput());
//  defconv->SetSpacing(timg->GetSpacing());
    defconv->Update();
    return defconv->GetOutput();

    }
  else
    {

    defreader->Update();
    return defreader->GetOutput();

    }

}


