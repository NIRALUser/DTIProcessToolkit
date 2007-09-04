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
