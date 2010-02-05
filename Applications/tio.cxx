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
    std::cerr << "Usage: " << argv[0] << " doffile imagefile outputfile [inputtype]" << std::endl;
    std::cerr << "with inputtype = " << std::endl;
    std::cerr << "0: old rview" << std::endl;
    std::cerr << "1: new rview (output txt file of dof2mat)" << std::endl;
    return EXIT_FAILURE;
    }

  bool rview_old = true;
  if(argc == 5)
  {
    if(atoi(argv[4]) == 1)
      rview_old = false;
    else if(atoi(argv[4]) == 0)
      rview_old = true;
    else
    {
      std::cerr << "inputtype has to be set to 0 (old rview) or 1 (new rview)" << std::endl;
      return EXIT_FAILURE;
    }
  }

  const std::string dofname(argv[1]);

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(argv[2]);
  reader->UpdateOutputInformation();

  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  ImageType::PointType origin = reader->GetOutput()->GetOrigin();


  typedef itk::AffineTransform<Precision, 3> AffineTransformType;

  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer twriter = TransformWriterType::New();

  //Handle the old version of rview (the one where the dof files were in ASCII)
  if(rview_old)
  {
    RViewTransform<Precision> dof = readDOFFile<Precision>(dofname);
    AffineTransformType::Pointer aff = createITKAffine(dof,
						       size,
						       spacing,
						       origin);
    twriter->AddTransform(aff);    
    twriter->SetFileName(argv[3]);
    twriter->Update();
  }
  //Handle the new version of rview (where the dof files are binary files and have to be converted into txt files with dof2mat)
  else
  {
    newRViewTransform<Precision> dof = readDOF2MATFile<Precision>(dofname);
    typedef itk::AffineTransform<Precision, 3> AffineTransformType;
    AffineTransformType::Pointer aff = createnewITKAffine(dof,
							  size,
							  spacing,
							  origin);
    twriter->AddTransform(aff);
    twriter->SetFileName(argv[3]);
    twriter->Update();
  } 
}
