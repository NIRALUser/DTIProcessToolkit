/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2007-09-05 19:35:36 $
  Version:   $Revision: 1.2 $
  Author:    Casey Goodlett (gcasey@sci.utah.edu)
             Pierre Fillard (Pierre.Fillard@sophia.inria.fr)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>
#include <vector>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageFileReader.h>
#include <itkDiffusionTensor3D.h>

#include "itkLogEuclideanTensorImageFilter.h"

int main (int argc, char* argv[])
{
  if(argc != 4)
    {
    std::cout << "Usage: " << argv[0]  << " file_in file_out sigma" << std::endl;
    return EXIT_FAILURE;
    }

  const char* file_in = argv[1];
  const char* file_out = argv[2];
  const double sigma = atof(argv[3]);

  typedef double                               ScalarType;
  typedef itk::DiffusionTensor3D<ScalarType>   TensorType;
  typedef itk::Image<TensorType, 3>            TensorImageType;
  typedef itk::Vector<double, 6>               VectorType;
  typedef itk::Image<VectorType, 3>            VectorImageType;

  typedef itk::Image<ScalarType, 3>            ImageType;

  std::cout << "Reading: " << file_in;
  std::cout << std::endl;
  
  // read in the tensor field
  typedef itk::ImageFileReader<TensorImageType> TensorReaderType;
  TensorReaderType::Pointer myIO = TensorReaderType::New();
  myIO->SetFileName (file_in);
  try
  {
    myIO->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e;
    return EXIT_FAILURE;
  }

  std::cout << "Done." << std::endl;

  TensorImageType::Pointer myTensorImage = myIO->GetOutput();


  std::cout << "Loging...";
  std::cout << std::endl;
  
  // log the tensor field
  typedef itk::LogEuclideanTensorImageFilter<double> LogFilterType;
  LogFilterType::Pointer myLoger = LogFilterType::New();
  myLoger->SetInput (myTensorImage);
  try
  {
    myLoger->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e;
    return EXIT_FAILURE;
  }

  std::cout << "Done." << std::endl;
  
  VectorImageType::Pointer myVectorImage = myLoger->GetOutput();
  
  std::cout << "Converting...";
  std::cout << std::endl;
  
  std::vector<ImageType::Pointer> ImageVector;

  // convert the tensor image into 6 scalar images
  for( unsigned int i=0 ;i<6; i++)
  {
    ImageType::Pointer image = ImageType::New();

    TensorImageType::RegionType region = myVectorImage->GetLargestPossibleRegion();

    image->SetRegions (region);
    image->SetSpacing (myVectorImage->GetSpacing());
    image->SetOrigin (myVectorImage->GetOrigin());
    image->SetDirection (myVectorImage->GetDirection());

    image->Allocate();

    itk::ImageRegionConstIteratorWithIndex<VectorImageType>  itIn (myVectorImage, myVectorImage->GetLargestPossibleRegion());
    itk::ImageRegionIteratorWithIndex<ImageType> itOut(image, image->GetLargestPossibleRegion());

    while( !itIn.IsAtEnd() )
    {

      VectorType T = itIn.Get();

      itOut.Set (T[i]);
      
      ++itIn;
      ++itOut;
    }

    
    ImageVector.push_back (image);
    
  }

  std::cout << "Done." << std::endl;

  std::cout << "Hessianing...";
  std::cout << std::endl;
  // now filters by computing the Hessian of each image
  typedef itk::HessianRecursiveGaussianImageFilter<ImageType>
    HessianFilterType;
  typedef HessianFilterType::OutputImageType HessianImageType;
  typedef HessianImageType::PixelType        HessianPixelType;
  
  std::vector<HessianImageType::Pointer> HessianImageVector;
  
  for( unsigned int i=0; i<6; i++)
  {
    HessianFilterType::Pointer myHessian = HessianFilterType::New();
    myHessian->SetInput (ImageVector[i]);
    myHessian->SetSigma (sigma);

    try
    {
      myHessian->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return -1;
    }

    HessianImageVector.push_back (myHessian->GetOutput());
  }

  std::cout << "Done." << std::endl;

  ImageVector.clear();

  std::cout << "Combining...";
  std::cout << std::endl;
  // finally: combines everything and extract the max eigenvalue  
  
  ImageType::Pointer myFinalImage  = ImageType::New();
  //ImageType::Pointer myFinalImage2 = ImageType::New();
  //ImageType::Pointer myFinalImage3 = ImageType::New();
  TensorImageType::RegionType region = myVectorImage->GetLargestPossibleRegion();
  
  myFinalImage->SetRegions (region);
  myFinalImage->SetSpacing (myVectorImage->GetSpacing());
  myFinalImage->SetOrigin (myVectorImage->GetOrigin());
  myFinalImage->SetDirection (myVectorImage->GetDirection());
  myFinalImage->Allocate();
  /*
  myFinalImage2->SetRegions (region);
  myFinalImage2->SetSpacing (myVectorImage->GetSpacing());
  myFinalImage2->SetOrigin (myVectorImage->GetOrigin());
  myFinalImage2->SetDirection (myVectorImage->GetDirection());
  myFinalImage2->Allocate();

  myFinalImage3->SetRegions (region);
  myFinalImage3->SetSpacing (myVectorImage->GetSpacing());
  myFinalImage3->SetOrigin (myVectorImage->GetOrigin());
  myFinalImage3->SetDirection (myVectorImage->GetDirection());
  myFinalImage3->Allocate();
  */

  typedef itk::ImageRegionConstIteratorWithIndex<HessianImageType> HessianImageIteratorType;
  
  std::vector<HessianImageIteratorType> HessianIteratorList;

  for( unsigned int i=0; i<6; i++)
  {
    HessianImageIteratorType it (HessianImageVector[i], HessianImageVector[i]->GetLargestPossibleRegion());
    HessianIteratorList.push_back (it);
  }

  itk::ImageRegionIteratorWithIndex<ImageType> itOut (myFinalImage, myFinalImage->GetLargestPossibleRegion());
  //itk::ImageRegionIteratorWithIndex<ImageType> itOut2 (myFinalImage2, myFinalImage2->GetLargestPossibleRegion());
  //itk::ImageRegionIteratorWithIndex<ImageType> itOut3 (myFinalImage3, myFinalImage3->GetLargestPossibleRegion());

  while( !itOut.IsAtEnd() )
  {

    TensorType SuperT;
    TensorType SuperT2;

    vnl_matrix<double> HH (3,18,0.0);
    
    for( unsigned int i=0; i<6; i++)
    {
      HessianPixelType H = HessianIteratorList[i].Get();

      // create a tensor struct for norm
//       TensorType T;
//       T[0] = H[0];
//       T[1] = H[1];
//       T[2] = H[3];
//       T[3] = H[2];
//       T[4] = H[4];
//       T[5] = H[5];

      HH (0,i*3+0) = H[0];
      HH (1,i*3+0) = H[1];
      HH (2,i*3+0) = H[3];
      HH (0,i*3+1) = H[1];
      HH (1,i*3+1) = H[2];
      HH (2,i*3+1) = H[4];
      HH (0,i*3+2) = H[3];
      HH (1,i*3+2) = H[4];
      HH (2,i*3+2) = H[5];
     
      /*
      double max = -999999999;
      for( unsigned int j=0; j<6; j++)
      {
        if ( T[j] > max )
        {
          max = T[j];
        }
      }
      
      SuperT[i]  = T.GetNorm();
      SuperT2[i] = T.GetEigenvalue (2);
      */
      //SuperT2[i] = max;
    }


    vnl_svd<double> MySVD (HH);

    //itOut.Set (SuperT.GetEigenvalue (2));
    //itOut2.Set (SuperT2.GetEigenvalue (2));
    itOut.Set (MySVD.W (0));


    ++itOut;
    //++itOut2;
    //++itOut3;
    for( unsigned int i=0; i<6; i++)
    {
      ++(HessianIteratorList[i]);
    }
    
  }

  std::cout << "Done." << std::endl;


  std::cout << "Writing...";
  std::cout << std::endl;



  std::string sfile_out = file_out;
  std::string::size_type pos = sfile_out.rfind (".gipl");
  if( pos != std::string::npos )
  {
    typedef itk::Image<unsigned short, 3> LightImageType;
    
    itk::RescaleIntensityImageFilter<ImageType, LightImageType>::Pointer rescaler=
      itk::RescaleIntensityImageFilter<ImageType, LightImageType>::New();
    
    rescaler->SetOutputMinimum ( 0 );
    rescaler->SetOutputMaximum ( 32767 );
    rescaler->SetInput ( myFinalImage );
    try
    {
      rescaler->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return -1;
    }

    itk::ImageFileWriter<LightImageType>::Pointer writer =
      itk::ImageFileWriter<LightImageType>::New();
    writer->SetFileName (file_out);
    writer->SetInput (rescaler->GetOutput());
    try
    {
      writer->Write();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return -1;
    }
    
  }
  else
  {
    // write the result:
    itk::ImageFileWriter<ImageType>::Pointer writer = itk::ImageFileWriter<ImageType>::New();
    writer->SetFileName (file_out);
    writer->SetInput (myFinalImage);
    try
    {
      writer->Write();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return -1;
    }

  }

  
  /*
  itk::ImageFileWriter<ImageType>::Pointer writer2 = itk::ImageFileWriter<ImageType>::New();
  std::string out2 = s_fileout + "_maxev.nrrd";
  writer2->SetFileName (out2.c_str());
  writer2->SetInput (myFinalImage2);
  try
  {
    writer2->Write();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e;
    return -1;
  }


  itk::ImageFileWriter<ImageType>::Pointer writer3 = itk::ImageFileWriter<ImageType>::New();
  std::string out3 = s_fileout + "_tensorev.nrrd";
  writer3->SetFileName (out3.c_str());
  writer3->SetInput (myFinalImage3);
  try
  {
    writer3->Write();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e;
    return -1;
  }
  */
  
  std::cout << "Done." << std::endl;
  
  
  return 0;
}

