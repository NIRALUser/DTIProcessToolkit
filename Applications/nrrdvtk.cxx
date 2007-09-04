#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>

#include "itkTensorFileReader.h"
#include "itkTensorFileWriter.h"

int main(int argc, char* argv[])
{
  typedef double RealType;
  typedef itk::DiffusionTensor3D<RealType> TensorType;
  typedef itk::Image<TensorType,3> TensorImageType;

  if(argc != 3)
    {
    std::cerr << "Usage: " << argv[0] << "input-file output-file" << std::endl;
    std::cerr << "This program converts between nrrd tensor files and vtk tensor files" << std::endl;
    return EXIT_FAILURE;
    }

  bool nrrd2vtk = true;
 
  const std::string file1(argv[1]);
  const std::string file2(argv[2]);

  if(file1.find(".vtk") != std::string::npos)
    {
    nrrd2vtk = false;
    }
  
  // vtk to nrrd
  if(!nrrd2vtk)
    {
    if(file2.find(".nhdr") == std::string::npos &&
       file2.find(".nrrd") == std::string::npos)
      {
      std::cerr << "Not really vtk to nrrd" << std::endl;
      }
    // do stuff
    typedef itk::TensorFileReader<TensorImageType> VTKTensorReaderType;
    VTKTensorReaderType::Pointer tread = VTKTensorReaderType::New();
    tread->SetFileName(file1.c_str());
    tread->Update();

    typedef itk::ImageFileWriter<TensorImageType> NRRDTensorWriterType;
    NRRDTensorWriterType::Pointer twrit = NRRDTensorWriterType::New();
    twrit->SetInput(tread->GetOutput());
    twrit->SetFileName(file2.c_str());
    twrit->Update();

    }
  // nrrd to vtk
  else
    {
    if(file2.find(".vtk") == std::string::npos)
      {
      std::cerr << "Not really nrrd to vtk" << std::endl;
      }
    typedef itk::ImageFileReader<TensorImageType> NRRDTensorReaderType;
    NRRDTensorReaderType::Pointer tread = NRRDTensorReaderType::New();
    tread->SetFileName(file1.c_str());
    tread->Update();

    typedef itk::TensorFileWriter<RealType> VTKTensorWriterType;
    VTKTensorWriterType::Pointer twrit = VTKTensorWriterType::New();
    twrit->SetInput(tread->GetOutput());
    twrit->SetFileName(file2.c_str());
    twrit->Update();

    }

  return EXIT_SUCCESS;
}