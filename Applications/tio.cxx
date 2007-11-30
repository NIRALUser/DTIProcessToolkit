#include <itkAffineTransform.h>
#include <itkTransformFileWriter.h>
#include <itkImageFileReader.h>
#include <string>

#include "transforms.h"

int main(int argc, char* argv[])
{
  typedef float Precision;
  typedef itk::Image<unsigned short,3 > ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  
  if(argc < 4)
    {
    std::cerr << "Usage: " << argv[0] << " doffile imagefile outputfile" << std::endl;
    return EXIT_FAILURE;
    }

  const std::string dofname(argv[1]);

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(argv[2]);
  reader->UpdateOutputInformation();

  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  ImageType::PointType origin = reader->GetOutput()->GetOrigin();

  RViewTransform<Precision> dof = readDOFFile<Precision>(dofname);

  typedef itk::AffineTransform<Precision, 3> AffineTransformType;
  AffineTransformType::Pointer aff = createITKAffine(dof,
                                                     size,
                                                     spacing,
                                                     origin);
  

  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer twriter = TransformWriterType::New();
  
  twriter->AddTransform(aff);
  twriter->SetFileName(argv[3]);
  
  twriter->Update();
}
