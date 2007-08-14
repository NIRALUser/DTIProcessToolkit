#include <iostream>
#include <fstream>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkChangeLabelImageFilter.h>

int main(int argc, char* argv[])
{
  if(argc != 4)
    {
    std::cerr << "Usage <infile> <outfile> <textfile>" << std::endl;
    return EXIT_FAILURE;
    }
  
  const char* infile = argv[1];
  const char* outfile = argv[2];
  const char* tablename = argv[3];

  typedef unsigned short PixelType;
  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageFileWriterType;

  typedef itk::ChangeLabelImageFilter<ImageType, ImageType> ChangeLabelType;

  ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
  reader->SetFileName(infile);

  ChangeLabelType::Pointer labelchanger = ChangeLabelType::New();
  labelchanger->SetInput(reader->GetOutput());

  std::ifstream tablereader(tablename);
  do
    {
    PixelType in, out;
    tablereader >> in >> out;
    std::cout << in << " to " << out << std::endl;
    labelchanger->SetChange(in, out);
    tablereader >> std::ws;
    }
  while(!tablereader.eof());

  ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
  writer->SetUseCompression(true);
  writer->SetInput(labelchanger->GetOutput());
  writer->SetFileName(outfile);
  writer->Update();

  return EXIT_SUCCESS;
}
