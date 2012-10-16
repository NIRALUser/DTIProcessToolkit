#include <iostream>
#include "itkRicianNoiseLevelDeterminer.h"
#include "itkImageFileReader.h"

const unsigned int DIM = 3;

typedef unsigned short DWIPixelType;
typedef float RealType;

typedef itk::Image<DWIPixelType, DIM> IntImageType;
typedef itk::ImageFileReader<IntImageType> FileReaderType;

typedef itk::RicianNoiseLevelDeterminer< IntImageType, RealType > RicianNoiseLevelDeterminerType;

int main()
{

  FileReaderType::Pointer fileReader = FileReaderType::New();
  fileReader->SetFileName( "baselineTst.nhdr" );

  try
    {
    std::cout << "Reading input file ... ";
    fileReader->Update();
    std::cout << "done." << std::endl;
    }
  catch (itk::ExceptionObject e)
    {
    std::cerr << " exception caught !" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }


  // now determine the Rician noise level

  RicianNoiseLevelDeterminerType::Pointer ricianNoiseLevelDeterminer = RicianNoiseLevelDeterminerType::New();
  
  ricianNoiseLevelDeterminer->SetInput( fileReader->GetOutput() );
  ricianNoiseLevelDeterminer->Compute();

  std::cout << "The Rician noise level is = " << ricianNoiseLevelDeterminer->GetOutput() << std::endl;

  return 0;
}
