#include <itkAffineTransform.h>
#include <itkTransformFileWriter.h>
#include <itkImageFileReader.h>
#include <string>
#include <itkPoint.h>
#include <itkVector.h>
#include "transforms.h"
#include <sstream>

void GetImageCenter( itk::Image< unsigned short , 3 >::Pointer image ,
                     itk::Image< unsigned short , 3 >::PointType &center
                   )
{
  typedef itk::Image< unsigned short , 3 > ImageType ;
  //Get lower corner position
  ImageType::PointType origin ;
  origin = image->GetOrigin() ;
  //Get higher corner position
  ImageType::SizeType size ;
  size = image->GetLargestPossibleRegion().GetSize() ;
  ImageType::IndexType index ;
  for( int i = 0 ; i < 3 ; i++ )
  {
    index[ i ] = size[ i ] - 1 ;
  }
  ImageType::PointType corner ;
  image->TransformIndexToPhysicalPoint( index , corner ) ;
  //Compute image center position
  for( int i = 0 ; i < 3 ; i++ )
  {
    center[ i ] = ( corner[ i ] + origin[ i ] ) / 2.0 ;
  }
}

template < class Precision >
int ComputeTransform(  const std::string doffile ,
                       std::string sourceFileName ,
                       const std::string targetFileName ,
                       const std::string outputFileName ,
                       bool rview_old
                     )
{
  typedef itk::Image< unsigned short , 3 > ImageType ;
  typedef itk::ImageFileReader< ImageType > ImageReaderType ;
  ImageReaderType::Pointer reader = ImageReaderType::New() ;
  reader->SetFileName( targetFileName.c_str() ) ;
  reader->UpdateOutputInformation() ;
  //Get the target image information to set the center of transform
  //to the center of the target image
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing() ;
  ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize() ;
  ImageType::PointType origin = reader->GetOutput()->GetOrigin() ;

  typedef itk::AffineTransform< Precision , 3 > AffineTransformType ;
  typedef itk::TransformFileWriter TransformWriterType ;
  TransformWriterType::Pointer twriter = TransformWriterType::New() ;
  typename AffineTransformType::Pointer aff ;
  //Handle the old version of rview (the one where the dof files were in ASCII)
  if( rview_old )
  {
    RViewTransform< Precision > dof = readDOFFile< Precision >( doffile ) ;
    aff = createITKAffine( dof ,
                           size ,
                           spacing ,
                           origin
                         ) ;
  }
  //Handle the new version of rview (where the dof files are binary
  //files and have to be converted into txt files with dof2mat)
  else
  {
    newRViewTransform< Precision > dof = readDOF2MATFile< Precision >( doffile ) ;
    typedef itk::AffineTransform< Precision , 3 > AffineTransformType ;
    aff = createnewITKAffine( dof,
							   size ,
							   spacing ,
							   origin
                                                         ) ;
  }
 //If the source image is given, add a translation that moves
 //the center of the target image to the center of the source image
 if( sourceFileName.compare( "" ) )
 {
   ImageReaderType::Pointer sourceReader = ImageReaderType::New() ;
   sourceReader->SetFileName( sourceFileName.c_str() ) ;
   sourceReader->UpdateOutputInformation() ;
   ImageType::PointType targetCenter ;
   GetImageCenter( reader->GetOutput() , targetCenter ) ;
   ImageType::PointType sourceCenter ;
   GetImageCenter( sourceReader->GetOutput() , sourceCenter ) ;
   itk::Vector< float , 3 > translation ;
   translation = sourceCenter - targetCenter ;
   translation += aff->GetTranslation() ;
   aff->SetTranslation( translation ) ;
 }
 twriter->AddTransform( aff ) ;    
 twriter->SetFileName( outputFileName ) ;
 try
 {
   twriter->Update() ;
 }
 catch(...)
 {
   return EXIT_FAILURE ;
 }
 return EXIT_SUCCESS ;
}



int main(int argc, char* argv[])
{
  typedef float Precision;
  typedef double doublePrecision ;
  std::string sourceFileName ;
  //Check if source is specified
  bool sourceSet = false ;
  int numberOfArgs = argc ;
  int outFilePos = 3 ;
  if( argc >= 6 && !strcmp( argv[ 3 ] , "-s") )
    {
      sourceSet = true ;
      numberOfArgs -= 2 ;
      sourceFileName.assign( argv[ 4 ] ) ;
      outFilePos = 5 ;
    }
  //check if old or new rview
  bool rview_old = true;
  int transformTypePos = outFilePos;
  if( argc > outFilePos + 1 )
  {
    std::istringstream inputType ;
    inputType.str( argv[ outFilePos + 1 ] ) ; 
    inputType>> rview_old ;
    if( inputType.fail() )
    {
      transformTypePos = outFilePos + 1 ;
    }
    else
    {
      numberOfArgs-- ;
      transformTypePos = outFilePos + 2 ;
    }
  }
  //Check if float or double
  bool doubleSet = false ;
  if( argc > transformTypePos && !strcmp( argv[ transformTypePos ] , "-d" ) )
  {
    doubleSet = true ;
    numberOfArgs-- ;
  }

  if( numberOfArgs != 4 )
    {
    std::cerr << "Usage: " << argv[0] << " doffile targetimagefile [-s sourceImageFileName] outputfile [inputtype] [-d (outputTransform as doubles instead of floats) ]" << std::endl;
    std::cerr << "with inputtype = " << std::endl;
    std::cerr << "0: new rview (output txt file of dof2mat)" << std::endl;//this was inverted
    std::cerr << "1: old rview" << std::endl;
    return EXIT_FAILURE;
    }

  const std::string doffile( argv[ 1 ] ) ;
  const std::string targetFileName( argv[ 2 ] ) ;
  const std::string outputFileName( argv[ outFilePos ] ) ;

  if( doubleSet )
  {
    return ComputeTransform< doublePrecision > ( doffile ,
                                                 sourceFileName ,
                                                 targetFileName ,
                                                 outputFileName ,
                                                 rview_old
                                               ) ;
  }
  else
  {
    return ComputeTransform< Precision > ( doffile ,
                                           sourceFileName ,
                                           targetFileName ,
                                           outputFileName ,
                                           rview_old
                                         ) ;
  }
}
