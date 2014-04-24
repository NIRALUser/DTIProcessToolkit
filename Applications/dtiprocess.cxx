/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/21 21:11:49 $
  Version:   $Revision: 1.7 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// STL includes
#include <string>
#include <iostream>

// ITK includes
// datastructures
#include <itkMetaDataObject.h>
#include <itkVersion.h>
// Filters
#include <itkCastImageFilter.h>
// IO
#include <itkImageFileReader.h>

// dtiprocess headers
#include "itkVectorMaskImageFilter.h"
#include "transforms.h"
#include "tensoroperations.h"
#include "imageio.h"
#include "deformationfieldio.h"

#include "dtiprocessCLP.h"

// Bad global variables.  TODO: remove these
bool VERBOSE = false;

int main(int argc, char* argv[])
{
  PARSE_ARGS;

  typedef itk::DiffusionTensor3D<float> TensorFloatPixelType;
  typedef itk::Image<TensorFloatPixelType, DIM> TensorFloatImageType;
  typedef itk::ImageFileWriter<TensorFloatImageType> TensorFileWriterType;
  typedef itk::CastImageFilter< TensorImageType, TensorFloatImageType > CastDTIFilterType ;
  // If the value scale is true (default) we scal FA and MD values to
  // integer ranges.
//   {
//     scale = false;
//   }
  bool scale = !scalarFloat; // I don't understand this logic but...

//   if(vm.count("verbose"))
//   {
//     VERBOSE = true;
//   }
  VERBOSE = verbose;
  // Read tensor image
  typedef itk::ImageFileReader<TensorImageType> FileReaderType;
  FileReaderType::Pointer dtireader = FileReaderType::New();
  if( dtiImage == "" )
    {
    std::cerr << "Missing DTI Image filename" << std::endl;
    exit(1);
    }
  dtireader->SetFileName(dtiImage.c_str() );
  try
    {
    dtireader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  // Parse gradient directions from image header as specified by the
  // namic conventions defined at http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
  GradientListType::Pointer gradientContainer = GradientListType::New();

  itk::MetaDataDictionary & dict = dtireader->GetOutput()->GetMetaDataDictionary();

  std::vector<std::string> keys = dict.GetKeys();
  for( std::vector<std::string>::const_iterator it = keys.begin();
       it != keys.end(); ++it )
    {
    if( it->find("DWMRI_gradient") != std::string::npos )
      {
      std::string value;

      itk::ExposeMetaData<std::string>(dict, *it, value);
      std::istringstream iss(value);
      GradientType       g;
      iss >> g[0] >> g[1] >> g[2];

      unsigned int ind;
      std::string  temp = it->substr(it->find_last_of('_') + 1);
      ind = atoi(temp.c_str() );

      gradientContainer->InsertElement(ind, g);
      }
    }
  for( std::vector<std::string>::const_iterator it = keys.begin();
       it != keys.end(); ++it )
    {
    if( it->find("DWMRI_NEX") != std::string::npos )
      {
      std::string numrepstr;

      itk::ExposeMetaData<std::string>(dict, *it, numrepstr);
      unsigned int numreps = atoi(numrepstr.c_str() );

      std::string  indtorepstr = it->substr(it->find_last_of('_') + 1);
      unsigned int indtorep =  atoi(indtorepstr.c_str() );

      GradientType g = gradientContainer->GetElement(indtorep);
      for( unsigned int i = indtorep + 1; i < indtorep + numreps; i++ )
        {
        gradientContainer->InsertElement(i, g);
        }
      }

    }

  // Debugging
  if( VERBOSE )
    {
    std::cout << "Interpolation type: "    // vm["interpolation"].as<InterpolationType>() << std::endl;
              << interpolation << std::endl;
    std::cout << "reorientation type: "    // vm["reorientation"].as<TensorReorientationType>() << std::endl;
              << reorientation << std::endl;
    }

  TensorImageType::Pointer tensors = dtireader->GetOutput();
  //  if(vm.count("mask"))
  if( mask != "" )
    {
    typedef itk::ImageFileReader<LabelImageType> MaskFileReaderType;
    MaskFileReaderType::Pointer maskreader = MaskFileReaderType::New();
    //    maskreader->SetFileName(vm["mask"].as<std::string>());
    maskreader->SetFileName(mask.c_str() );
    typedef itk::VectorMaskImageFilter<TensorImageType, LabelImageType, TensorImageType> MaskFilterType;
    MaskFilterType::Pointer _mask = MaskFilterType::New();
    _mask->SetInput1(tensors);
    _mask->SetInput2(maskreader->GetOutput() );

    try
      {
      _mask->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }

    tensors = _mask->GetOutput();
    // If the outmask option is specified, the masked tensor field is saved
    if( outmask != "" )
      {
      if( !doubleDTI )
        {
        CastDTIFilterType::Pointer castFilter = CastDTIFilterType::New() ;
        castFilter->SetInput( tensors ) ;
        TensorFileWriterType::Pointer dtiwriter = TensorFileWriterType::New();
        dtiwriter->SetFileName(outmask.c_str());
        dtiwriter->SetInput(castFilter->GetOutput());
        dtiwriter->SetUseCompression(true);
        try
          {
          dtiwriter->Update();
          }
          catch (itk::ExceptionObject & e)
          {
            std::cerr << e <<std::endl;
            return EXIT_FAILURE;
          }
        }
      else
        {
        typedef itk::ImageFileWriter<TensorImageType> FileWriterType;
        FileWriterType::Pointer dtiwriter = FileWriterType::New();
        dtiwriter->SetFileName(outmask.c_str());
        dtiwriter->SetUseCompression(true);
        dtiwriter->SetInput(tensors);
        try
          {
          dtiwriter->Update();
          }
        catch (itk::ExceptionObject & e)
          {
          std::cerr << e <<std::endl;
          return EXIT_FAILURE;
          }
        }
      }
    }

  // sigma set in PARSE_ARGS
  // double sigma = vm["sigma"].as<double>();

  // Compute FA image
  //  if(vm.count("fa-output"))
  if( faOutput != "" )
    {
    if( scale )
      {
      writeImage(faOutput,
                 createFA<unsigned short>(tensors) );
      }
    else
      {
      writeImage(faOutput,
                 createFA<double>(tensors) );
      }
    }

  //  if(vm.count("fa-gradient-output"))
  if( faGradientOutput != "" )
    {
    writeImage(faGradientOutput,
               createFAGradient(tensors, sigma) );
    }

  //  if(vm.count("fa-gradmag-output"))
  if( faGradientMagOutput != "" )
    {
    writeImage(faGradientMagOutput,
               createFAGradMag(tensors, sigma) );
    }

  if( colorFAOutput != "" )
    {
    writeImage(colorFAOutput,
               createColorFA(tensors) );
    }

  if( principalEigenvectorOutput != "" )
//     vm.count("closest-dotproduct-output"))
    {
    writeImage(principalEigenvectorOutput,
               createPrincipalEigenvector(tensors) );
    }

  //  if(vm.count("md-output"))
  if( mdOutput != "" )
    {
    if( scale )
      {
      writeImage(mdOutput,
                 createMD<unsigned short>(tensors) );
      }
    else
      {
      writeImage(mdOutput,
                 createMD<double>(tensors) );
      }
    }

  if( lambda1Output != "" )
    {
    if( scale )
      {
      writeImage(lambda1Output,
                 createLambda<unsigned short>(tensors, Lambda1) );
      }
    else
      {
      writeImage(lambda1Output,
                 createLambda<double>(tensors, Lambda1) );
      }
    }

  if( lambda2Output != "" )
    {
    if( scale )
      {
      writeImage(lambda2Output,
                 createLambda<unsigned short>(tensors, Lambda2) );
      }
    else
      {
      writeImage(lambda2Output,
                 createLambda<double>(tensors, Lambda2) );
      }
    }

  if( lambda3Output != "" )
    {
    if( scale )
      {
      writeImage(lambda3Output,
                 createLambda<unsigned short>(tensors, Lambda3) );
      }
    else
      {
      writeImage(lambda3Output,
                 createLambda<double>(tensors, Lambda3) );
      }
    }

  if( RDOutput != "" )
    {
    if( scale )
      {
      writeImage(RDOutput,
                 createRD<unsigned short>(tensors) );
      }
    else
      {
      writeImage(RDOutput,
                 createRD<double>(tensors) );
      }
    }

  if( frobeniusNormOutput != "" )
    {
    if( scale )
      {
      writeImage(frobeniusNormOutput,
                 createFro<unsigned short>(tensors) );
      }
    else
      {
      writeImage(frobeniusNormOutput,
                 createFro<double>(tensors) );
      }
    }

  if( negativeEigenvectorOutput != "" )
    {
    writeImage(negativeEigenvectorOutput,
               createNegativeEigenValueLabel(tensors) );
    }

  if( rotOutput != "" )
    {
    TensorImageType::Pointer tensorImage ;
    //If the input affine file is a dof file from rview
    if(dofFile != "")
      {
      tensorImage = createROT(tensors,dofFile,0) ;
      }
    //If the input affine file is a new dof file (output of dof2mat)
    else if(newdof_file != "")
      {
      tensorImage = createROT(tensors, newdof_file, 1);
      }
    //If the input affine file is an itk compatible file
    else if(affineitk_file != "")
      {
      tensorImage = createROT(tensors, affineitk_file, 2);
      }
    else
      {
      std::cerr << "Tensor rotation requested, but dof/newdof/affineitk file not specified" << std::endl;
      return EXIT_FAILURE;
      }
    if( !doubleDTI )
      {
      CastDTIFilterType::Pointer castFilter = CastDTIFilterType::New() ;
      castFilter->SetInput( tensorImage ) ;
      castFilter->Update() ;
      TensorFloatImageType::Pointer tensorFloat = castFilter->GetOutput() ;
      writeImage( rotOutput , tensorFloat ) ;
      }
    else
      {
      writeImage( rotOutput , tensorImage ) ;
      }
    }
  if( deformationOutput != "" )
    {
    if( forwardTransformation == "" )
      {
      std::cerr << "Deformation field info not fully specified" << std::endl;
      return EXIT_FAILURE;
      }
    DeformationImageType::Pointer forward;

    DeformationFieldType dftype = Displacement;
    if( hField )
      {
      dftype = HField;
      }

    forward = readDeformationField(forwardTransformation, dftype);
    TensorImageType::Pointer tensorImage ;
    tensorImage = createWarp(tensors,
               forward,
               //vm["reorientation"].as<TensorReorientationType>(),
               (reorientation == "fs" ? FiniteStrain :
               PreservationPrincipalDirection),
               //vm["interpolation"].as<InterpolationType>()));
               (interpolation == "linear" ? Linear :
               (interpolation == "nearestneightbor" ? NearestNeighbor :
               Cubic)));
    if( !doubleDTI )
      {
      CastDTIFilterType::Pointer castFilter = CastDTIFilterType::New() ;
      castFilter->SetInput( tensorImage ) ;
      castFilter->Update() ;
      TensorFloatImageType::Pointer tensorFloat = castFilter->GetOutput() ;
      writeImage( deformationOutput , tensorFloat ) ;
      }
    else
      {
      writeImage( deformationOutput , tensorImage ) ;
      }
    }

#if 0 //
  if( vm.count("stats") )
    {
    if( !vm.count("mask") )
      {
      std::cerr << "WARNING: No mask specified.  Computing whole brain statistics" << std::endl;
      }
    if( VERBOSE )
      {
      std::cout << "Computing tensor stats" << std::endl;
      }

//    TensorRegionStatistics tstat = computeTensorStatistics(tensors);

    }
#endif
  return EXIT_SUCCESS;
}
