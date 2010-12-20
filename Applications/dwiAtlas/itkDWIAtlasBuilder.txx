/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkDWIAtlasBuilder.txx,v $
Language:  C++
Date:      $Date: 2007-01-26 09:45:07 $
Version:   $Revision: 1.8 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDWIAtlasBuilder_txx
#define __itkDWIAtlasBuilder_txx

#include "itkDWIAtlasBuilder.h"
#include <itkNumericTraits.h>
#include "itkProgressReporter.h"

namespace itk
{

/** Constructor */
template <class MyRealType, class DWIPixelType>
DWIAtlasBuilder<MyRealType, DWIPixelType>
::DWIAtlasBuilder()
{
  this->SetNumberOfRequiredInputs(1);
  m_Size.Fill(0);
  m_Spacing.Fill(1.0);
  m_Origin.Fill(0.0);

  this->nrOfDatasets = 0;

  this->m_Verbose = true;
  this->m_IsHField = true;

  this->m_GradientVectorFile = "None";
  this->m_MaskImageFileName = "None";

  this->m_SHOrder = -1;
  this->m_Lambda = 0;
  this->m_JustDoResampling = false;
  this->m_NoLogFit = false;

  this->m_DoWeightedLS = false;

  this->m_InterpolationType = NEAREST_NEIGHBOR;
  this->m_AveragingType = GEOMETRIC;

  this->m_ScalingFactor = 1.0;
  this->m_LogMinArgumentValue = 1;

  this->m_HuberC = 2.0;
  this->m_RiceSigma = -1;  // needs to be data-specific, otherwise meaningless
  this->m_NrOfWLSIterations = 3;

  this->m_UsedOrder = 0;
  this->m_NumTerms = 0;
  this->m_NrOfBaselines = 0;
  this->m_numnewgvectors = 0;

  this->dwiFiles.clear();
  this->deformationFiles.clear();

  this->dwireader = NULL;
  this->deformation = NULL;
  this->gradientContainers = NULL;
  this->jacobian = NULL;

  this->m_NumberOfThreads = -1;

  // TODO: fix this

  this->m_ConsoleProgressCommandPointer = ConsoleProgressCommand::New();
  /*this->AddObserver( StartEvent(), m_ConsoleProgressCommandPointer );
  this->AddObserver( EndEvent(), m_ConsoleProgressCommandPointer );
  this->AddObserver( ProgressEvent(), m_ConsoleProgressCommandPointer );*/

}

/** Destructor */
template <class MyRealType, class DWIPixelType>
DWIAtlasBuilder<MyRealType, DWIPixelType>
::~DWIAtlasBuilder()
{
  delete [] dwireader;
  delete [] deformation;
  delete [] gradientContainers;
  delete [] jacobian;
}

template <class MyRealType, class DWIPixelType>
void
DWIAtlasBuilder<MyRealType, DWIPixelType>
::SetNrOfThreads( int iNrOfThreads )
{
  if ( iNrOfThreads>0 )
    {
    m_NumberOfThreads = iNrOfThreads;
    }
}

template <class MyRealType, class DWIPixelType>
unsigned int
DWIAtlasBuilder<MyRealType, DWIPixelType>
::GetNrOfThreads( void )
{
  return this->GetNumberOfThreads();
}

template <class MyRealType, class DWIPixelType>
void
DWIAtlasBuilder<MyRealType, DWIPixelType>
::SetInterpolationType( std::string interpolationType )
{
  if ( interpolationType.compare( "nearestNeighbor" )==0 )
    {
    m_InterpolationType = NEAREST_NEIGHBOR;
    }
  else if ( interpolationType.compare( "linearlyInterpolated" )==0 )
    {
    m_InterpolationType = LINEAR;
    }
  else if ( interpolationType.compare( "weightedInterpolation" )==0 )
    {
    m_InterpolationType = USE_ALL_WITH_WEIGHTING;
    }
  else
    {
    m_InterpolationType = NEAREST_NEIGHBOR;
    std::cout << "WARNING: unknown interpolation type = " << interpolationType << std::endl;
    std::cout << "Defaulting to " << this->GetInterpolationType() << std::endl;
    }
}

template <class MyRealType, class DWIPixelType>
std::string
DWIAtlasBuilder<MyRealType, DWIPixelType>
::GetInterpolationType( void )
{
  switch ( m_InterpolationType )
    {
    case NEAREST_NEIGHBOR:
      return "nearestNeighbor";
      break;
    case LINEAR:
      return "linearlyInterpolated";
      break;
    case USE_ALL_WITH_WEIGHTING:
      return "weightedInterpolation";
      break;
    default:
      std::cout << "WARNING: Unkown interpolation type " << m_InterpolationType << std::endl;
      std::cout << "You should never see this message. Resetting to nearest neighbor." << std::endl;
      m_InterpolationType = NEAREST_NEIGHBOR;
      return "unknownInterpolationType::WARNING";
    }
}

template <class MyRealType, class DWIPixelType>
void
DWIAtlasBuilder<MyRealType, DWIPixelType>
::SetAveragingType( std::string averagingType )
{
  if ( averagingType.compare( "algebraic" )==0 )
    {
    m_AveragingType = ALGEBRAIC;
    }
  else if ( averagingType.compare( "geometric" )==0 )
    {
    m_AveragingType = GEOMETRIC;
    }
  else
    {
    m_AveragingType = GEOMETRIC;
    std::cout << "WARNING: unknown averaging type = " << averagingType << std::endl;
    std::cout << "Defaulting to " << this->GetAveragingType() << std::endl;
    }
}

template <class MyRealType, class DWIPixelType>
std::string
DWIAtlasBuilder<MyRealType, DWIPixelType>
::GetAveragingType( void )
{
  switch ( m_AveragingType )
    {
    case ALGEBRAIC:
      return "algebraic";
      break;
    case GEOMETRIC:
      return "geometric";
      break;
    default:
      std::cout << "WARNING: Unkown averaging type " << m_AveragingType << std::endl;
      std::cout << "You should never see this message. Resetting to geometric." << std::endl;
      m_AveragingType = GEOMETRIC;
      return "unknownInterpolationType::WARNING";
    }
}
 
template <class MyRealType, class DWIPixelType>
void
DWIAtlasBuilder<MyRealType, DWIPixelType>
::LoadDataAndInitialize()
{

  // initialize global structures for faster spherical harmonics
  // computations; important (!), will crash otherwise
  
   sh::initSH();

   // set the case list (input parameter) and prepare everything to
   // load the data

   aux::getFiles( *(this->GetInput()->Get() ), dwiFiles, deformationFiles, m_JustDoResampling ); 
   this->nrOfDatasets = dwiFiles.size();

   // read in all the DWI data

   std::vector<std::string>::iterator iterDwiFiles;
   std::vector<std::string>::iterator iterDeformationFiles;

   dwireader = new typename FileReaderType::Pointer[ nrOfDatasets ];

   unsigned int cCount = 0;
   for ( iterDwiFiles=dwiFiles.begin(); iterDwiFiles!=dwiFiles.end(); iterDwiFiles ++ )
     {
     dwireader[ cCount ] = FileReaderType::New();
     if ( m_Verbose )
       std::cout << "Setting dwi filename to " << *iterDwiFiles << std::endl;
     dwireader[ cCount ]->SetFileName( *iterDwiFiles );
     cCount++;
     }

   // read in all the deformation fields

   deformation = new DeformationImageType::Pointer[ nrOfDatasets ];

   cCount = 0;
   
   if ( !m_JustDoResampling )
     {
     for ( iterDeformationFiles=deformationFiles.begin(); iterDeformationFiles!=deformationFiles.end(); iterDeformationFiles++ )
       { 
       if ( m_IsHField )
	 deformation[ cCount ] = readDeformationField( *iterDeformationFiles, HField);
       else
	 deformation[ cCount ] = readDeformationField( *iterDeformationFiles, Displacement );
       cCount++;
       }
     }
   
   // declare the gradient containers that will hold the gradient
   // information for all of the images (chances are they will all be
   // the same in practice, but let's keep this general for now)

   gradientContainers = new typename DiffusionEstimationFilterType::GradientDirectionContainerType::Pointer[ nrOfDatasets ];
   
   m_NrOfBaselines = 0;
   
   std::vector<unsigned int> vecBaselineIndices;

   for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
     {
     gradientContainers[iI] = DiffusionEstimationFilterType::GradientDirectionContainerType::New();
     itk::MetaDataDictionary & dict = dwireader[iI]->GetOutput()->GetMetaDataDictionary();
     if ( m_Verbose )
       {
       std::cout << "Parsing the meta dictionary of " << dwiFiles[iI] << std::endl;
       }
     dwireader[iI]->Update(); // otherwise it does not know what the
     // meta data dictionary contains
     // convert the direction matrix to the currently set real type
     vnl_matrix<MyRealType> currentDirectionMatrix(3,3,0);
     vnl_matrix<double> dwiReaderDirectionMatrix = dwireader[iI]->GetOutput()->GetDirection().GetVnlMatrix();
     for ( int iX=0; iX<3; iX++ )
       {
       for ( int iY=0; iY<3; iY++ )
	 {
	 currentDirectionMatrix(iX,iY) = (MyRealType)dwiReaderDirectionMatrix(iX,iY);
	 }
       }
     
     getGradients( *gradientContainers[iI], dict, currentDirectionMatrix, m_NrOfBaselines, vecBaselineIndices, m_Verbose );

     if ( m_DoWeightedLS && m_RiceSigma<=0 && (vecBaselineIndices.size()>0) )
       {
       // need to compute an estimation for the noise level in the
       // data

       typename RicianNoiseLevelDeterminerType::Pointer ricianNoiseLevelDeterminer = RicianNoiseLevelDeterminerType::New();
  
       typename ExtractInputVolumeFilterType::Pointer eviFilter = ExtractInputVolumeFilterType::New();
  
       if ( m_Verbose ) std::cout << "Estimating noise level ...";

       eviFilter->SetInput( this->dwireader[iI]->GetOutput() );
       eviFilter->SetVolumeNr( vecBaselineIndices[0] );
       eviFilter->Update();

       ricianNoiseLevelDeterminer->SetInput( eviFilter->GetOutput() );
       ricianNoiseLevelDeterminer->Compute();

       if ( m_Verbose ) std::cout << "done." << std::endl;

       m_RiceSigma = ricianNoiseLevelDeterminer->GetOutput();

       if ( m_Verbose )
	 {
	 std::cout << "Noise level = " << m_RiceSigma << std::endl;
	 }
       }

     }


   // get the desired gradients or clone them

   if ( m_GradientVectorFile.compare( "None" )==0 )
     {
     // no gradient file specified, so just use the same gradients as
     // in the first input volume
     
     // go through all of them (of the first volume) and put them into
     // the desired gradient file; baselines will be passed through
     // averaged
     
     // first find out how many gradients there are
     
     // the number of the DWIs are all we have in the gradient
     // containter - the number of baselines
     
     if ( m_Verbose )
       std::cout << "Cloning direction information from input" << std::endl;
     
     unsigned int iNrNonBaselineDWIs = gradientContainers[0]->Size() - m_NrOfBaselines;
     m_DesiredGradients.set_size( iNrNonBaselineDWIs, 3 );
     
     unsigned int iCurrentGradient = 0;
     
     for ( unsigned int iI=0; iI<gradientContainers[0]->Size(); iI++ )
       {
       vnl_vector_fixed<MyRealType, 3> tgrad = gradientContainers[0]->ElementAt(iI);
       if ( tgrad.magnitude()!=0 )
	 {
	 m_DesiredGradients.set_row( iCurrentGradient++, tgrad );
	 }
       }
     
     }
   else
     {
     // read them from the file
     unsigned int iNrNonBaselineDWIs = aux::getNumOutputGradients( m_GradientVectorFile );
     m_DesiredGradients.set_size( iNrNonBaselineDWIs, 3 );
     aux::parseGradientFile<MyRealType>( m_GradientVectorFile, m_DesiredGradients );
     }
   
  // construct the spherical harmonics matrix for the desired output gradients
   
  //convert newgrads to spherical coordinates
   m_numnewgvectors = m_DesiredGradients.rows();
   vnl_matrix<MyRealType> newgvectorssph(m_numnewgvectors, 2);
   sh::cart2sph<MyRealType>(m_DesiredGradients, newgvectorssph);
   
   // number of gradients is the number of datasets multiplied by the
   // number of non-baseline DWIs per dataset. Note, it is not
   // advisable to use the highest possible order for the spherical
   // harmonics expansion, since the gradient directions will in most
   // cases be strongly aligned (unless there are substantial
   // registration warps occuring
   
   unsigned int totalNumberOfGradientsInAtlasSpace = (gradientContainers[0]->Size()-m_NrOfBaselines)*nrOfDatasets; 

   m_NumTerms = sh::getNumTermsAndCheckOrder( m_SHOrder, totalNumberOfGradientsInAtlasSpace, m_UsedOrder );
   
   std::cout << "Used order = " << m_UsedOrder << std::endl;
   
   //construct SHBnew
   m_sh_basis_mat_new.set_size(m_numnewgvectors, m_NumTerms);
   m_sh_basis_mat_new.fill( 0 );
   sh::generateSHBasisMatrix<MyRealType>(m_NumTerms, m_numnewgvectors, newgvectorssph, m_sh_basis_mat_new);

   // and compute the Jacobians

   jacobian = new typename MyJacobianFilterType::Pointer[ nrOfDatasets ];

   if ( !m_JustDoResampling )
     {
     for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
       {
       jacobian[iI] = MyJacobianFilterType::New();
       // Compute Jacobian of transform
       jacobian[iI]->SetUseImageSpacingOn();
       jacobian[iI]->SetInput(deformation[iI]);
       jacobian[iI]->Update();
       }
     }

   if ( m_Verbose )
     std::cout << "Number of components per pixel (input) = " << dwireader[0]->GetOutput()->GetNumberOfComponentsPerPixel() << std::endl;

   // if we have a mask image, load it, otherwise create one of the
   // appropriate size (this avoids code replication)

   if ( m_MaskImageFileName.compare("None")!=0 )
     {
     // we have one, load it

     maskReader = ScalarFileReaderType::New();
     maskReader->SetFileName( m_MaskImageFileName );

     if ( m_Verbose )
       std::cout << "Reading mask image from " << m_MaskImageFileName << " ... ";

     maskReader->Update();

     if ( m_Verbose )
       std::cout << "done." << std::endl;

     MaskImage = maskReader->GetOutput();

     }
   else
     {
     // we do not have one create one
     MaskImage = ScalarImageType::New();
     MaskImage->SetOrigin( dwireader[0]->GetOutput()->GetOrigin() );
     MaskImage->SetSpacing( dwireader[0]->GetOutput()->GetSpacing() );
     MaskImage->SetRegions( dwireader[0]->GetOutput()->GetLargestPossibleRegion() );
     MaskImage->Allocate();

     // initialize it with all ones

     typedef itk::ImageRegionIterator<ScalarImageType> ScalarImageIteratorType;
     ScalarImageIteratorType maskIterator( MaskImage, dwireader[0]->GetOutput()->GetLargestPossibleRegion() );

     maskIterator.GoToBegin();

     while ( !maskIterator.IsAtEnd() )
       {
       maskIterator.Set( 1 );
       ++maskIterator;
       }

     }
   
}

template <class MyRealType, class DWIPixelType>
vnl_matrix<MyRealType>
DWIAtlasBuilder<MyRealType, DWIPixelType>
::rotateGradients( typename GradientDirectionContainerType::Pointer gradientContainer, itk::Matrix<MyRealType, 3, 3> jacobian, std::vector<bool>& isBaseline, unsigned int &iNrOfBaselines, bool bypassRotation )
{
  const unsigned int numgrads = gradientContainer->Size(); 
  vnl_matrix_fixed<MyRealType, 3, 3> iden;
  iden.set_identity();
  vnl_matrix_fixed<MyRealType, 3, 3> localt = vnl_inverse(jacobian.GetVnlMatrix() + iden);
  

  // use polar decompostion to get the rotation matrix out of localt
  typedef vnl_matrix<MyRealType> VnlMatrixType;

  VnlMatrixType PQ = localt;
  VnlMatrixType NQ = localt;
  VnlMatrixType PQNQDiff;

  const unsigned int maximumIterations = 100;

  for(unsigned int ni = 0; ni < maximumIterations; ni++ )
     {
     // Average current Qi with its inverse transpose
     NQ = ( PQ + vnl_inverse_transpose( PQ ) ) / 2.0;
     PQNQDiff = NQ - PQ;
     if( PQNQDiff.frobenius_norm() < 1e-7 )
       {
       if ( m_Verbose ){
        std::cout << "Polar decomposition used "
                  << ni << " iterations " << std::endl;
        break;
        }
       }
     else
       {
       PQ = NQ;
       }
     }

  VnlMatrixType QMatrix;

  QMatrix = NQ;
  if ( m_Verbose )
    {
      std::cout << "Initial Matrix = " << std::endl << localt << std::endl;
      std::cout << "Q Matrix = " << std::endl << QMatrix << std::endl;
    }

  iNrOfBaselines = 0;
  isBaseline.clear();

  vnl_matrix<MyRealType> rotatedGradients( numgrads ,3 );

  for (unsigned int iM = 0; iM < numgrads; iM++)
  {
    vnl_vector_fixed<MyRealType, 3> tgrad = QMatrix * gradientContainer->ElementAt(iM);
    if(tgrad.magnitude())
      {
      tgrad /= tgrad.magnitude();
      isBaseline.push_back( false );
      }
    else
      {
      iNrOfBaselines++;
      isBaseline.push_back( true );
      }
    
    /*if ( m_Verbose )
      std::cout << "Found " << iNrOfBaselines << " baseline(s)." << std::endl;*/

    rotatedGradients.set_row( iM, tgrad );

  }
  return rotatedGradients;

}

template <class MyRealType, class DWIPixelType>
void 
DWIAtlasBuilder<MyRealType, DWIPixelType>
::getGradients(  typename DiffusionEstimationFilterType::GradientDirectionContainerType & gradientContainer, itk::MetaDataDictionary & dict,  vnl_matrix<MyRealType> imgf, unsigned int& iNrOfBaselines, std::vector<unsigned int> &vecBaselineIndices, bool m_Verbose )
{

  vecBaselineIndices.clear();

  // Read dwi meta-data

  // read into bvalue the DWMRI_b-value that is the b-value of
  // the experiment
  MyRealType bvalue = 0;
  bool readbvalue = false;
  
  iNrOfBaselines = 0;

  // read into gradientContainer the gradients

  //if ( m_Verbose )
  //  dict.Print( std::cout );

  vnl_matrix<MyRealType> transform(3,3);
  transform.set_identity();

  // Apply measurement frame if it exists

  // TODO: MAYBE PUT BACK IN?

  std::cout << "WARNING: IGNORING ANY MEASUREMENT FRAME INFORMATION" << std::endl;
  

  /*if(dict.HasKey(NRRD_MEASUREMENT_KEY))
  {
    if(m_Verbose)
      std::cout << "Reorienting gradient directions to image coordinate frame" << std::endl;

    // measurement frame
    vnl_matrix<MyRealType> mf(3,3);
    // imaging frame

    //std::vector<std::vector<MyRealType> > nrrdmf;
    //itk::ExposeMetaData<std::vector<std::vector<MyRealType> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);

    std::vector<std::vector<double> > nrrdmf;
    itk::ExposeMetaData<std::vector<std::vector<double> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);

    //if(m_Verbose)
    //{
    //  std::cout << "Image frame: " << std::endl;
    //  std::cout << imgf << std::endl;
    //}

    for(unsigned int i = 0; i < 3; ++i)
    {
      for(unsigned int j = 0; j < 3; ++j)
      {
      mf(i,j) = (MyRealType)nrrdmf[j][i];

      nrrdmf[j][i] = (double)imgf(i,j);
      }
    }

    //if(m_Verbose)
    //{
    //  std::cout << "Measurement frame: " << std::endl;
    //  std::cout << mf << std::endl;
    //}

    //itk::EncapsulateMetaData<std::vector<std::vector<MyRealType> >
    //>(dict,NRRD_MEASUREMENT_KEY,nrrdmf);

    itk::EncapsulateMetaData<std::vector<std::vector<double> > >(dict,NRRD_MEASUREMENT_KEY,nrrdmf);

    // prevent slicer error

    transform = vnl_svd<MyRealType>(imgf).inverse()*mf;

    //if(m_Verbose)
    //{
    //  std::cout << "Transform: " << std::endl;
    //  std::cout << transform << std::endl;
    //}

  }*/

  if(dict.HasKey("modality"))
  {
    itk::EncapsulateMetaData<std::string>(dict,"modality","DWMRI");
    if ( m_Verbose ) std::cout << "modality = DWMRI"<< std::endl;
  }

  std::vector<std::string> keys = dict.GetKeys();

  for(std::vector<std::string>::const_iterator it = keys.begin();
      it != keys.end(); ++it)
    {
    std::string value;
    if( it->find("DWMRI_b-value") != std::string::npos)
      {
      std::string t;
      itk::ExposeMetaData<std::string>(dict, *it, t);
      readbvalue = true;
      bvalue = atof(t.c_str());
      }
    else if( it->find("DWMRI_gradient") != std::string::npos)
      {
      std::string value;

      itk::ExposeMetaData<std::string>(dict, *it, value);
      std::istringstream iss(value);
      MyGradientType g;
      iss >> g[0] >> g[1] >> g[2];
      
      g = transform * g;
      
      unsigned int ind;
      std::string temp = it->substr(it->find_last_of('_')+1);
      ind = atoi(temp.c_str());
      
      if ( g.magnitude()==0 )
	{
	vecBaselineIndices.push_back( ind );
	iNrOfBaselines++;
	}
      
      gradientContainer.InsertElement(ind,g);
      }
    else if( it->find("DWMRI_NEX") != std::string::npos)
      {
      std::string numrepstr;
      
      itk::ExposeMetaData<std::string>(dict, *it, numrepstr);
      unsigned int numreps = atoi(value.c_str());
      
      std::string indtorepstr = it->substr(it->find_last_of('_')+1);
      unsigned int indtorep =  atoi(indtorepstr.c_str());
      
      MyGradientType g = gradientContainer.GetElement(indtorep);
      
      for(unsigned int i = indtorep+1; i < indtorep+numreps-1; i++)
        gradientContainer.InsertElement(i,g);

      if ( g.magnitude()==0 )
	iNrOfBaselines+=numreps;
    }

  }

  if(!readbvalue)
  {
    std::cerr << "BValue not specified in header file" << std::endl;
    //return EXIT_FAILURE; // TODO: fix this, so it can really
    //abort/return something meaningful
  }

  if(m_Verbose)
    std::cout << "BValue: " << bvalue << std::endl;

}

template <class MyRealType, class DWIPixelType>
void
DWIAtlasBuilder<MyRealType, DWIPixelType>
::computeAveragedBaselines( std::vector<MyRealType>& averagedBaselines, STransformedGradientInformationType& allBaselines, unsigned int iNrOfBaselines, unsigned int interpolationType, unsigned int averagingType, unsigned int &nrOfBaselineOutliers )
{
  // computes an average baseline based on the extracted baseline data
  // default is to create as many baselines as the original dataset
  // (though this could be easily adapted in a later version of the
  // program)
  //
  // iNrOfBaselines is the number of baselines of the original image

  nrOfBaselineOutliers = 0;

  averagedBaselines.clear();
  averagedBaselines.resize( iNrOfBaselines, 0 );

  switch ( interpolationType )
    {
    case NEAREST_NEIGHBOR:
    case LINEAR:
      // the same code applies for nearest neighbor and linear
    {
    unsigned int iNrOfMeasurementsForAveraging = allBaselines.dwiVals.size()/iNrOfBaselines; 
    
    if ( averagingType==ALGEBRAIC )
      {
      // loop over each complete baseline sets
      for ( unsigned int iI=0; iI<iNrOfMeasurementsForAveraging; iI++ )
	{
	// loop over one set and add its contribution
	for ( unsigned int iJ=0; iJ<iNrOfBaselines; iJ++ )
	  {
	  averagedBaselines[iJ] += allBaselines.dwiVals[iI*iNrOfBaselines + iJ];
	  }
	}
      // take the average, by dividing by the numer of elements used
      for ( unsigned int iI=0; iI<iNrOfBaselines; iI++ )
	{
	averagedBaselines[iI]/= iNrOfMeasurementsForAveraging;
	}
      }
    else if ( averagingType==GEOMETRIC )
      {
      // do the same (as for algebraic) but with geometric averaging, i.e., sum the
      // logs then take the exponential
      
      if ( !m_DoWeightedLS ) //just compute the geometric average
	{

	// loop over each complete baseline sets
	for ( unsigned int iI=0; iI<iNrOfMeasurementsForAveraging; iI++ )
	  {
	  // loop over one set and add its log contribution
	  for ( unsigned int iJ=0; iJ<iNrOfBaselines; iJ++ )
	    {
	    averagedBaselines[iJ] += log( (MyRealType)std::max( allBaselines.dwiVals[iI*iNrOfBaselines + iJ], m_LogMinArgumentValue) );
	    }
	  }
	// take the average, by dividing by the numer of elements used
	for ( unsigned int iI=0; iI<iNrOfBaselines; iI++ )
	  {
	  averagedBaselines[iI]/= iNrOfMeasurementsForAveraging;
	  }
	
	// now take the exponential
	
	for ( unsigned int iI=0; iI<iNrOfBaselines; iI++ )
	  {
	  // this is now the geometric mean: exp(1/n\sum_{i=1}^n log S_i)
	  averagedBaselines[iI] = exp( averagedBaselines[iI] ); 
	  }
	}
      else // do weighted least squares, by iteration, supports robust function
	{
	
	// get initial value by computing an unweighted version as
	// above

	std::vector<MyRealType> huberWeights;
	huberWeights.resize( iNrOfMeasurementsForAveraging );

	std::vector<MyRealType> logY;
	logY.resize( iNrOfBaselines*iNrOfMeasurementsForAveraging );

	std::vector<MyRealType> logSE;
	logSE.resize( iNrOfBaselines );

	// first log transform all the baseline images

	for ( unsigned int iI=0; iI<iNrOfBaselines*iNrOfMeasurementsForAveraging; iI++ )
	  {
	  logY[iI] = log( (MyRealType)std::max( allBaselines.dwiVals[iI], m_LogMinArgumentValue) );
	  }

	// compute initial estimates for the log transformed
	// measurements SE, by not using any weighting

	for ( unsigned int iJ=0; iJ<iNrOfBaselines; iJ++ )
	  {
	  logSE[iJ] = 0;
	  for ( unsigned int iI=0; iI<iNrOfMeasurementsForAveraging; iI++ )
	    {
	    logSE[iJ] += logY[iI*iNrOfBaselines + iJ];
	    }
	  logSE[iJ] /= iNrOfMeasurementsForAveraging;
	  }


	// given these intitial conditions, we can iterate the
	// weighted least squares solution

	// now iterate to obtain updated versions for all of these
	// do it averaged volume by averaged volume
	
	bool isOutlier;

	for ( unsigned int iJ=0; iJ<iNrOfBaselines; iJ++ ) 
	  { // loop over all the volumes
	  unsigned int iNrOfOutliersOfCurrentAveraging;
	  for ( unsigned int iW=0; iW<m_NrOfWLSIterations; iW++ )
	    {
	    iNrOfOutliersOfCurrentAveraging = 0;
	    // compute the Huber coefficients for the current estimate
	    // and the current set of measurements
	    MyRealType currentSE = exp( logSE[iJ] );
	    for ( unsigned int iI=0; iI<iNrOfMeasurementsForAveraging; iI++ )
	      {
	      MyRealType scaledResidual = (logSE[iJ]-logY[iI*iNrOfBaselines + iJ])*currentSE/m_RiceSigma;
	      huberWeights[iI] = sh::huberWeightFcn( scaledResidual, m_HuberC, isOutlier );
	      huberWeights[iI] *= currentSE*currentSE/(m_RiceSigma*m_RiceSigma);

	      if ( isOutlier ) iNrOfOutliersOfCurrentAveraging++;

	      }

	    // now compute new estimate of the logs based on this

	    logSE[iJ] = 0;
	    MyRealType huberWeightSum = 0;
	    for ( unsigned int iI=0; iI<iNrOfMeasurementsForAveraging; iI++ )
	      {
	      logSE[iJ] += huberWeights[iI]*logY[iI*iNrOfBaselines + iJ];
	      huberWeightSum += huberWeights[iI];
	      }
	    logSE[iJ]/=huberWeightSum;

	    } // end for iW

	  nrOfBaselineOutliers += iNrOfOutliersOfCurrentAveraging;

	  }

	// now we have estimates and can get back to the real value by
	// computing the exponent

	for ( unsigned int iJ=0; iJ<iNrOfBaselines; iJ++ ) 
	  {
	  averagedBaselines[iJ] = exp( logSE[iJ] ); 
	  }


	} // else do weighted least squares
      
      } // geometric averaging

    } // FIXME: Not sure why these extrac curly braces are needed. Is
      // there a code issue somewhere?
    break;
    
    case USE_ALL_WITH_WEIGHTING:

      // TODO!!

      /*unsigned int iNrOfMeasurementsForAveraging = allBaselines.dwiVals.size()/(iNrOfBaselines*NCONTROLPOINTS);

      MyRealType totalWeight = 0;
      
      if ( averagingType==ALGEBRAIC )
	{
	
	}
      else if ( averagingType==GEOMETRIC )
	{
	}*/

      break;
      // only the weighted case is different
      // we do the same looping, but will have more elements based on
      // the number of control points
    }


}

template <class MyRealType, class DWIPixelType>
void 
DWIAtlasBuilder<MyRealType, DWIPixelType>
::extractDesiredBaselinesAndDWIs( STransformedGradientInformationType& allBaselines, STransformedGradientInformationType& allDWIs, const STransformedGradientInformationType* transformedInformation, unsigned int nrOfDataSets , unsigned int interpolationType, unsigned int averagingType )
{
  // implements different methods of information extraction
  // ignores all measurements which would require the use of invalid
  // points

  // first determine the number of usable datapoints (i.e., the ones
  // that did not result in an out-of-domain interpolation)

  //bool bSUPERVERBOSE = true;

  unsigned int iNrOfUsableDataPoints = 0;
  unsigned int iTotalNumberOfUsableBaselinesPerVoxel = 0;
  unsigned int iTotalNumberOfUsableDWIsPerVoxel = 0;

  bool foundAGoodVoxel = false;

  for ( unsigned int iI=0; iI<nrOfDataSets; iI++ )
    {
    if ( transformedInformation[iI].goodVoxel ) 
      {
      // TODO: Check that this works properly if we do not do interpolation
      iNrOfUsableDataPoints++;
      iTotalNumberOfUsableBaselinesPerVoxel += transformedInformation[iI].iNrOfBaselinesPerVolume;
      iTotalNumberOfUsableDWIsPerVoxel += transformedInformation[iI].NDWI-transformedInformation[iI].iNrOfBaselinesPerVolume;
      foundAGoodVoxel = true;
      }
    }

  /*if ( bSUPERVERBOSE )
    {
    if ( foundAGoodVoxel )
      std::cout << "good voxel found." << std::endl;
    }*/

  unsigned int iNumberOfStoredDWI = 0;
  unsigned int iNumberOfStoredBaseline = 0;

  switch ( interpolationType )
    {
    case NEAREST_NEIGHBOR:
      // for each of the gradient sets, find the one with the highest
      // weight value and use it

      allBaselines.gradients.set_size( iTotalNumberOfUsableBaselinesPerVoxel, DIM );
      allBaselines.gradients.fill ( 0 );  // all zero, because these are the baselines
      allBaselines.dwiVals.resize( iTotalNumberOfUsableBaselinesPerVoxel );
      allBaselines.isBaseline.resize( iTotalNumberOfUsableBaselinesPerVoxel, true );
      allBaselines.interpolationWeights.resize( iTotalNumberOfUsableBaselinesPerVoxel, 1 ); // uniform
      allBaselines.NDWI = iTotalNumberOfUsableBaselinesPerVoxel;
      allBaselines.iNrOfBaselinesPerVolume = iTotalNumberOfUsableBaselinesPerVoxel;
      if ( foundAGoodVoxel )
	allBaselines.goodVoxel = true;
      else
	allBaselines.goodVoxel = false;

      allDWIs.gradients.set_size( iTotalNumberOfUsableDWIsPerVoxel, DIM );
      allDWIs.dwiVals.resize( iTotalNumberOfUsableDWIsPerVoxel );
      allDWIs.isBaseline.resize( iTotalNumberOfUsableDWIsPerVoxel, false );
      allDWIs.interpolationWeights.resize( iTotalNumberOfUsableDWIsPerVoxel, 1 ); // uniform
      allDWIs.NDWI = iTotalNumberOfUsableDWIsPerVoxel;
      allDWIs.iNrOfBaselinesPerVolume = 0;
      if ( foundAGoodVoxel )
	allDWIs.goodVoxel = true;
      else
	allDWIs.goodVoxel = false;
      
      // now go through all of them and assign them based on their
      // highest weights

      for ( unsigned int iI=0; iI<nrOfDataSets; iI++ )
	{
	const unsigned int NDWI = transformedInformation[iI].NDWI;
	if ( transformedInformation[iI].goodVoxel )
	  {
	  typename std::vector<MyRealType>::const_iterator min_it = min_element( transformedInformation[iI].interpolationWeights.begin(), transformedInformation[iI].interpolationWeights.begin() );
	  // now get which index this corresponds to
	  unsigned int iIndex = min_it-transformedInformation[iI].interpolationWeights.begin();
	  
	  // now enter all the information at this index

	  for ( unsigned int iJ=0; iJ<NDWI; iJ++ )
	    {
	    if ( transformedInformation[iI].isBaseline[iJ] )
	      {
	      allBaselines.dwiVals[iNumberOfStoredBaseline] = transformedInformation[iI].dwiVals[iIndex*NDWI+iJ];
	      iNumberOfStoredBaseline++;
	      }
	    else
	      {
	      allDWIs.dwiVals[iNumberOfStoredDWI] = transformedInformation[iI].dwiVals[iIndex*NDWI+iJ];
	      allDWIs.gradients.set_row( iNumberOfStoredDWI, transformedInformation[iI].gradients.get_row( iIndex*NDWI+iJ ) );
	      iNumberOfStoredDWI++;
	      }
	    }


	  }
	}


      /*if ( iNumberOfStoredBaseline>0 )
	{
	std::cout << "number of stored baselines = " << iNumberOfStoredBaseline << std::endl;
	std::cout << "number of stored DWIs = " << iNumberOfStoredDWI << std::endl;
	}*/

      break;
    case LINEAR:
      std::cout << "Linear interpolation model not implemented! Abort." << std::endl;
      exit(-1);
      // does pretty much the same as for nearest neighbor but
      // averages pixels. This can be done with or without the rician
      // noise model. The number of resulting measurements will be the
      // same

      /*allBaselines.gradients.set_size( iTotalNumberOfUsableBaselinesPerVoxel, DIM );
      allBaselines.gradients.fill ( 0 );  // all zero, because these are the baselines
      allBaselines.dwiVals.resize( iTotalNumberOfUsableBaselinesPerVoxel );
      allBaselines.interpolationWeights.resize( iTotalNumberOfUsableBaselinesPerVoxel, 1 ); // uniform
      allBaselines.isBaseline.resize( iTotalNumberOfUsableBaselinesPerVoxel, true );
      allBaselines.NDWI = iTotalNumberOfUsableBaselinesPerVoxel;
      allBaselines.iNrOfBaselinesPerVolume = iTotalNumberOfUsableBaselinesPerVoxel;
      if ( foundAGoodVoxel )
	allBaselines.goodVoxel = true;
      else
	allBaselines.goodVoxel = false;

      allDWIs.gradients.set_size( iTotalNumberOfUsableDWIsPerVoxel, DIM );
      allDWIs.dwiVals.resize( iTotalNumberOfUsableDWIsPerVoxel );
      allDWIs.interpolationWeights.resize( iTotalNumberOfUsableDWIsPerVoxel, 1 ); // uniform
      allDWIs.isBaseline.resize( iTotalNumberOfUsableDWIsPerVoxel, false );
      allDWIs.NDWI = iTotalNumberOfUsableDWIsPerVoxel;
      allDWIs.iNrOfBaselinesPerVolume = 0;
      if ( foundAGoodVoxel )
	allDWIs.goodVoxel = true;
      else
	allDWIs.goodVoxel = false;
      
      // now go through all of them and weight them
      // either algebraically or geometrically

      // TODO: continue here

      for ( unsigned int iI=0; iI<nrOfDataSets; iI++ )
	{
	const unsigned int NDWI = transformedInformation[iI].NDWI;
	if ( transformedInformation[iI].goodVoxel )
	  {
	  
	  // now do the weighting
	  // set everything to zero first

	  iNumberOfStoredBaseline = 0;
	  iNumberOfStoredDWI = 0;

	  // initialize everything to zero first

	  for ( unsigned int iJ=0; iJ<NDWI; iJ++ )
	    {
	    if ( transformedInformation[iI].isBaseline[iJ] )
	      {
	      allBaselines.dwiVals[iNumberOfStoredBaseline] = 0;
	      iNumberOfStoredBaseline++;
	      }
	    else
	      {
	      allDWIs.dwiVals[iNumberOfStoredDWI] = 0;
	      iNumberOfStoredDWI++;
	      }
	    }
	  allDWIs.gradients.fill( 0 );  // intialize to zero gradients
	  
	  // now do the averaging

	  if ( averagingType==ALGEBRAIC ) 
	    {
	    // just add based on weighted sum
	    
	    iNumberOfStoredBaseline = 0;
	    iNumberOfStoredDWI = 0;
	    
	    for ( unsigned int iJ=0; iJ<NDWI; iJ++ )
	      {
	      if ( transformedInformation[iI].isBaseline[iJ] )
		{
		for ( unsigned int iK=0; iK<NCONTROLPOINTS; iK++ )
		  {
		  // FIMXE: cast is not appropriate
		  allBaselines.dwiVals[iNumberOfStoredBaseline] += transformedInformation[iI].interpolationWeights[iK]*transformedInformation[iI].dwiVals[iK*NDWI+iJ];
		  }
		iNumberOfStoredBaseline++;
		}
	      else
		{
		for ( unsigned int iK=0; iK<NCONTROLPOINTS; iK++ )
		  {
		  allDWIs.dwiVals[iNumberOfStoredDWI] += transformedInformation[iI].interpolationWeights[iK]*transformedInformation[iI].dwiVals[iK*NDWI+iJ];
		  // CHOICE: how do we average the gradient direction?
		  // Simply add, normalize afterwards
		  allDWIs.gradients.set_row( iNumberOfStoredDWI, allDWIs.gradients.get_row( iNumberOfStoredDWI ) + transformedInformation[iI].gradients.get_row( iK*NDWI+iJ ) );
		  }
		// now give the summed up gradient of norm one
		// TODO: add support for gradient directions which are
		// not unit length!
		
		MyRealType currentNorm = allDWIs.gradients.get_row( iNumberOfStoredDWI ).two_norm();
		allDWIs.gradients.scale_row( iNumberOfStoredDWI,  1.0/currentNorm );

		iNumberOfStoredDWI++;
		}
	      }
	    // there is no averaging necessary, because the weights
	    // are such that they sum to one! (CHECK THIS)
	    }
	  else if ( averagingType==GEOMETRIC )
	    {
	    // add geometrically, weight the logarithm then take the
	    // exponent

	    iNumberOfStoredBaseline = 0;
	    iNumberOfStoredDWI = 0;
	    
	    for ( unsigned int iJ=0; iJ<NDWI; iJ++ )
	      {
	      if ( transformedInformation[iI].isBaseline[iJ] )
		{
		for ( unsigned int iK=0; iK<NCONTROLPOINTS; iK++ )
		  {
		  allBaselines.dwiVals[iNumberOfStoredBaseline] += transformedInformation[iI].interpolationWeights[iK]*log( std::max( transformedInformation[iI].dwiVals[iK*NDWI+iJ], m_LogMinArgumentValue ) );
		  }
		iNumberOfStoredBaseline++;
		}
	      else
		{
		for ( unsigned int iK=0; iK<NCONTROLPOINTS; iK++ )
		  {
		  allDWIs.dwiVals[iNumberOfStoredDWI] += transformedInformation[iI].interpolationWeights[iK]*log( std::max( transformedInformation[iI].dwiVals[iK*NDWI+iJ], m_LogMinArgumentValue ) );
		  // CHOICE: how do we average the gradient direction?
		  // Simply add, normalize afterwards
		  allDWIs.gradients.set_row( iNumberOfStoredDWI, allDWIs.gradients.get_row( iNumberOfStoredDWI ) + transformedInformation[iI].gradients.get_row( iK*NDWI+iJ ) );
		  }
		// now give the summed up gradient of norm one
		// TODO: add support for gradient directions which are
		// not unit length!
		
		MyRealType currentNorm = allDWIs.gradients.get_row( iNumberOfStoredDWI ).two_norm();
		allDWIs.gradients.scale_row( iNumberOfStoredDWI,  1.0/currentNorm );

		iNumberOfStoredDWI++;
		}
	      }
	    // there is no averaging necessary, because the weights
	    // are such that they sum to one! (CHECK THIS)
	    // but we still need to take the exponent of all these
	    // values

	    iNumberOfStoredBaseline = 0;
	    iNumberOfStoredDWI = 0;

	    for ( unsigned int iJ=0; iJ<NDWI; iJ++ )
	      {
	      if ( transformedInformation[iI].isBaseline[iJ] )
		{
		allBaselines.dwiVals[iNumberOfStoredBaseline] = (DWIPixelType)round( exp( allBaselines.dwiVals[iNumberOfStoredBaseline] ) );
		iNumberOfStoredBaseline++;
		}
	      else
		{
		allDWIs.dwiVals[iNumberOfStoredDWI] = (DWIPixelType)round(exp( allDWIs.dwiVals[iNumberOfStoredDWI] ));
		iNumberOfStoredDWI++;
		}
	      }	    
	    }
	  }
	}*/
      break;
    case USE_ALL_WITH_WEIGHTING:
      // here we just add everything, split into baseline and dwis


/*      CONTINUE HERE: Check that the dimensions are correct!
	TODO: MAKE SURE TO INTRODUCE A CLEAR CONVENTION WHAT THE VARIABLES MEAN THAT IS CONSISTEN ACROSS THE DIFFERENT PROCESSING MODES

      allBaselines.gradients.set_size( iTotalNumberOfUsableBaselinesPerVoxel*NCONTROLPOINTS, DIM );
      allBaselines.gradients.fill ( 0 );  // all zero, because these are the baselines
      allBaselines.dwiVals.resize( iTotalNumberOfUsableBaselinesPerVoxel*NCONTROLPOINTS );
      allBaselines.interpolationWeights.resize( iTotalNumberOfUsableBaselinesPerVoxel*NCONTROLPOINTS, 1 ); // uniform
      allBaselines.isBaseline.resize( iTotalNumberOfUsableBaselinesPerVoxel*NCONTROLPOINTS, true );
      allBaselines.NDWI = iTotalNumberOfUsableBaselinesPerVoxel*NCONTROLPOINTS;
      allBaselines.iNrOfBaselinesPerVolume = iTotalNumberOfUsableBaselinesPerVoxel*NCONTROLPOINTS;
      if ( foundAGoodVoxel )
	allBaselines.goodVoxel = true;
      else
	allBaselines.goodVoxel = false;

      allDWIs.gradients.set_size( iTotalNumberOfUsableDWIsPerVoxel*NCONTROLPOINTS, DIM );
      allDWIs.dwiVals.resize( iTotalNumberOfUsableDWIsPerVoxel*NCONTROLPOINTS );
      allDWIs.interpolationWeights.resize( iTotalNumberOfUsableDWIsPerVoxel*NCONTROLPOINTS, 1 ); // uniform
      allDWIs.isBaseline.resize( iTotalNumberOfUsableDWIsPerVoxel*NCONTROLPOINTS, false );
      allDWIs.NDWI = iTotalNumberOfUsableDWIsPerVoxel*NROFCONTROLPOINTS;
      allDWIs.iNrOfBaselinesPerVolume = 0;
      if ( foundAGoodVoxel )
	allDWIs.goodVoxel = true;
      else
	allDWIs.goodVoxel = false;*/

      std::cout << "Weighting model not implemented! Abort." << std::endl;
      exit(-1);
      break;
    }
}




/** Set the Input caselist */
template <class MyRealType, class DWIPixelType>
void 
DWIAtlasBuilder<MyRealType, DWIPixelType>
::SetInput( const std::string* psCaseList )
{
  
  // strings are not dataobjects, so need to decorate it to push it down
  // the pipeline 
  typename InputStringObjectType::Pointer stringObject = 
    InputStringObjectType::New();
  stringObject->Set( const_cast< std::string*  >( psCaseList ) );
  
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0,  stringObject  );
}

template <class MyRealType, class DWIPixelType>
const typename DWIAtlasBuilder< MyRealType, DWIPixelType >::InputStringObjectType *
DWIAtlasBuilder<MyRealType, DWIPixelType>
::GetInput(void) 
{
  if (this->GetNumberOfInputs() < 1)
    {
    return 0;
    }
  return static_cast<const InputStringObjectType * >
    (this->ProcessObject::GetInput(0) );
}


template <class MyRealType, class DWIPixelType>
void 
DWIAtlasBuilder<MyRealType, DWIPixelType>
::SetSpacing(const double* spacing)
{
  SpacingType s(spacing);
  this->SetSpacing( s );
}


template <class MyRealType, class DWIPixelType>
void 
DWIAtlasBuilder<MyRealType, DWIPixelType>
::SetOrigin(const double* origin)
{
  PointType p(origin);
  SetOrigin(p);
}

template <class MyRealType, class DWIPixelType>
void
DWIAtlasBuilder<MyRealType, DWIPixelType>
::GenerateOutputInformation()
{
  // Get the input and output pointers 
  // Get from decorator
  OutputImageType     *outputImage    = this->GetOutput();

  // before we do this we need to know how big everything this, so
  // let's load the data now

  LoadDataAndInitialize();

  // Set output image params and Allocate image, do this by cloning
  // the information of the first dwi image (we assume all have the
  // same resolution -- otherwise there will be trouble)
  typename OutputImageType::RegionType region;
  region.SetSize( dwireader[0]->GetOutput()->GetLargestPossibleRegion().GetSize() );
  outputImage->SetRegions( region );

  // set spacing and origin

  outputImage->SetSpacing( dwireader[0]->GetOutput()->GetSpacing() );   
  outputImage->SetOrigin(  dwireader[0]->GetOutput()->GetOrigin() );   
		
									
  outputImage->SetNumberOfComponentsPerPixel( m_NrOfBaselines + m_numnewgvectors );
  outputImage->Allocate();

  // if we compute a weighted version, also allocate a volume
  // indicating the number of found outliers

  OutlierImage = ScalarImageType::New();
  OutlierImage->SetOrigin( dwireader[0]->GetOutput()->GetOrigin() );
  OutlierImage->SetSpacing( dwireader[0]->GetOutput()->GetSpacing() );
  OutlierImage->SetRegions( region );
  OutlierImage->Allocate();

  itk::MetaDataDictionary & dict = dwireader[0]->GetOutput()->GetMetaDataDictionary();
  
  // create a new metadata dictionary
  itk::MetaDataDictionary newDictionary;
  
  aux::ConstructOutputMetaDataDictionary<MyRealType>( newDictionary, dict, m_NrOfBaselines, m_DesiredGradients );
    
  outputImage->SetMetaDataDictionary( newDictionary );
}

//----------------------------------------------------------------------------

template< class MyRealType, class DWIPixelType >
void
DWIAtlasBuilder< MyRealType, DWIPixelType >
::AfterThreadedGenerateData( void )
{
  // so that the progress is on a new line
  std::cout << std::endl;
}

template< class MyRealType, class DWIPixelType >
void
DWIAtlasBuilder< MyRealType, DWIPixelType >
::BeforeThreadedGenerateData( void )
{
  if ( m_NumberOfThreads>0 )
    {
    std::cout << "Setting number of threads to " << m_NumberOfThreads << std::endl;
    this->SetNumberOfThreads( m_NumberOfThreads );
    std::cout << "Set number of threads  to = " << this->GetNumberOfThreads() << std::endl;
    }

  this->m_ConsoleProgressCommandPointer->SetMaxProgress( 100*this->GetNumberOfThreads() );

}

/** Update */
template< class MyRealType, class DWIPixelType >
void
DWIAtlasBuilder< MyRealType, DWIPixelType >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
  itkDebugMacro(<< "DWIAtlasBuilder::Update() called");

  // Get the output pointer (and allocates at the same time??)

  OutputImageType     *outputImage    = this->GetOutput();
 
  //ProgressReporter progress( this, threadId,
  //outputRegionForThread.GetNumberOfPixels() );
  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels(), 100 );

  /*std::cout << std::endl << "tid = " << threadId << " " << outputRegionForThread << " #of pixels = " << outputRegionForThread.GetNumberOfPixels() << std::endl;
  std::cout.flush();*/

  // setup the iterators over all images and all the deformation fields

  typedef itk::ImageRegionConstIterator<VectorImageType> VectorIterator;
  typedef itk::ImageRegionConstIteratorWithIndex<DeformationImageType> HFieldIterator;
  typedef itk::ImageRegionIterator<VectorImageType> dwisReconstructedIteratorType;
  typedef itk::ImageRegionIterator<ScalarImageType> OutlierImageIteratorType;
  typedef itk::ImageRegionConstIterator<ScalarImageType> MaskImageIteratorType;
  

  VectorIterator dwiits[ nrOfDatasets ];
  HFieldIterator hfieldits[ nrOfDatasets ];
  dwisReconstructedIteratorType recondwiit(outputImage, outputRegionForThread );

  OutlierImageIteratorType outlierit( OutlierImage, outputRegionForThread );
  MaskImageIteratorType maskit( MaskImage, outputRegionForThread );

  // iterate over all the deformation fields and the output image at
  // the same time

  recondwiit.GoToBegin();
  outlierit.GoToBegin();
  maskit.GoToBegin();

  for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
    {
    dwiits[iI] = VectorIterator(dwireader[iI]->GetOutput(), outputRegionForThread );
    dwiits[iI].GoToBegin();
    }
  
  if ( !m_JustDoResampling )
    {
    for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
      {
      hfieldits[iI] = HFieldIterator(deformation[iI], outputRegionForThread );
      hfieldits[iI].GoToBegin();
      }
    }

  bool noMoreElements = false;

  unsigned int iNrOfBaselines = 0; // TODO: check what we really want here

  // storage for the pseudo-inverse
  vnl_matrix<MyRealType> sh_basis_mat_pi; 

  // storage in case of weighted least squares

  vnl_matrix<MyRealType> B;
  vnl_diag_matrix<MyRealType> L;
  vnl_diag_matrix<MyRealType> W;

  unsigned long myCounter = 0;
  
  unsigned long nrOfPixels = outputRegionForThread.GetNumberOfPixels(); 

  /*if ( m_Verbose )
    {
    std::cout << "m_DoWeightedLS = " << m_DoWeightedLS << std::endl;
    std::cout << "m_HuberC = " << m_HuberC << std::endl;
    std::cout << "m_Lambda = " << m_Lambda << std::endl;
    std::cout << "m_RiceSigma = " << m_RiceSigma << std::endl;
    std::cout << "m_NrOfWLSIterations = " << m_NrOfWLSIterations << std::endl;
    }*/

  MyJacobianType j;
  j.SetIdentity();

  while ( !noMoreElements )
    {


    if ( myCounter%(nrOfPixels/20)==0 )
      {
      std::cout << "thread id = " << threadId << " " << 100*((double)myCounter)/nrOfPixels << "%" << std::endl;
      }
    myCounter++;

    // get values for all the DWI's at their respective locations
    // as given by the deformation field

    // loop over all the cases and get them one by one

    STransformedGradientInformationType transformedInformation[ nrOfDatasets ];

    bool isInMask = (maskit.Get()>0);

    if ( isInMask )
      {

      for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
	{
	
	itk::Index<3> pind;
	DeformationPixelType arggh;
	itk::ContinuousIndex<MyRealType, 3u> ci;
	
	transformedInformation[iI].dwiVals.clear();
	transformedInformation[iI].interpolationWeights.clear();
	transformedInformation[iI].isBaseline.clear();
	
	// TODO: Assume same spacing for now, adapt to different spacing
	// for images
	
	const typename VectorImageType::SpacingType spacing = dwireader[iI]->GetOutput()->GetSpacing();

	if ( m_JustDoResampling )
	  {
	  ci = maskit.GetIndex();
	  }
	else
	  {	// need to look at the deformation field
	  pind = hfieldits[iI].GetIndex();
	  arggh = hfieldits[iI].Get();
	  ci[0] = arggh[0] / spacing[0] + pind[0];
	  ci[1] = arggh[1] / spacing[1] + pind[1];
	  ci[2] = arggh[2] / spacing[2] + pind[2];
	  }
	
	// now get all the DWIs for this particular location and surrounding
	
	std::vector<MyRealType> weight(8, 1.0);
	MyRealType tweight = 0.0;
	
	// max, max, max
	// max, max, min
	// max, min, max
	// max, min, min
	// min, max, max
	// min, max, min
	// min, min, max
	// min, min, min
	transformedInformation[iI].goodVoxel = true;
	
	const unsigned int NDWI = dwiits[iI].Get().GetSize();
	transformedInformation[iI].NDWI = NDWI;
	
	transformedInformation[iI].gradients.set_size(NCONTROLPOINTS*NDWI,DIM);
	
	for(unsigned int iJ = 0; iJ < NCONTROLPOINTS; ++iJ)
	  {
	  
	  itk::Index<3u> cpi;
	  for(unsigned int mask = 0; mask < DIM; ++mask)
	    {
	    if(iJ & (1 << mask))
	      {
	      cpi[mask] = static_cast<long int>(floor(ci[mask]));
	      weight[iJ] *= 1 - (ci[mask] - floor(ci[mask]));
	      }
	    else
	      {
	      cpi[mask] = static_cast<long int>(ceil(ci[mask]));
	      weight[iJ] *= ci[mask] - floor(ci[mask]);
	      }
	    }
	  tweight += weight[iJ];
	  
	  if(!dwireader[iI]->GetOutput()->GetLargestPossibleRegion().IsInside(cpi))
	    {
	    transformedInformation[iI].goodVoxel = false;
	    break; // TODO: maybe take this out in case it interfers
	    // with the calculation of the baselines by dividing
	    // by the number of control points
	    }
	  
	  itk::VariableLengthVector<DWIPixelType> dwi = dwireader[iI]->GetOutput()->GetPixel(cpi);
	  
	  // store the dwi values
	  for (unsigned int iK=0; iK<NDWI; iK++)
	    {
	    transformedInformation[iI].dwiVals.push_back( dwi[iK] );
	    }
	  //std::copy(dwi.GetDataPointer(), dwi.GetDataPointer()+NDWI, transformedInformation[iI].dwiVals.begin() + iJ*NDWI );
	  
	  
	  vnl_matrix<MyRealType> rotatedGradientDirections;
	  
	  if ( m_JustDoResampling )
	    {
	    // do not apply the rotation, this can be used for debugging 
	    rotatedGradientDirections = rotateGradients( gradientContainers[iI], j, transformedInformation[iI].isBaseline, iNrOfBaselines, true );  
	    }
	  else
	    {
	    j = jacobian[iI]->GetOutput()->GetPixel(cpi);

	    rotatedGradientDirections = rotateGradients( gradientContainers[iI], j, transformedInformation[iI].isBaseline, iNrOfBaselines );
	    }
	  
	  transformedInformation[iI].iNrOfBaselinesPerVolume = iNrOfBaselines; 
	  
	  /*if ( m_Verbose )
	  std::cout << "Setting baselines to " << iNrOfBaselines << std::endl;*/
	  
	  // store the gradient directions into the matrix
	  
	  transformedInformation[iI].gradients.update( rotatedGradientDirections, iJ*NDWI,0 );
	  
	  } // end loop over control points

	// add the normalized weights
	
	for(unsigned int iJ = 0; iJ < NCONTROLPOINTS; ++iJ)
	  {
	  transformedInformation[iI].interpolationWeights.push_back( weight[iJ]/tweight );
	  }
	
	} /// end loop over datasets

      }

    // now we have all the information in the transformedInformation
    // data structure

    // extract the relevant subinformation, baselines and actual
    // gradients that are being used

    // here is also where interpolation happens or does not happen

    STransformedGradientInformationType allBaselines;
    STransformedGradientInformationType allDWIs;

    vnl_vector<MyRealType> s_values_new;
    vnl_vector<MyRealType> dwiVals;

    unsigned int nrOfDWIOutliers = 0;
    unsigned int nrOfBaselineOutliers = 0;
      

    if ( isInMask )
      {
      
      //unsigned int interpolationType = NEAREST_NEIGHBOR;
      //unsigned int interpolationType = LINEAR;
      //unsinged int interpolationType = USE_ALL_WITH_WEIGHTING;
      
      /*if ( m_Verbose )
      std::cout << "There are " << transformedInformation[0].NDWI << " dwis total and " << transformedInformation[0].iNrOfBaselinesPerVolume << " baselines per volume." << std::endl;*/
      
      extractDesiredBaselinesAndDWIs( allBaselines, allDWIs, transformedInformation, nrOfDatasets , m_InterpolationType, m_AveragingType );
      
      s_values_new.set_size( m_numnewgvectors );
      dwiVals.set_size( allDWIs.dwiVals.size() );

      // now let's compute the new approximation

      if ( m_NoLogFit )
	{
	for (unsigned int iV=0; iV<allDWIs.dwiVals.size(); iV++ ) dwiVals[iV] = allDWIs.dwiVals[iV];
	}
      else
	{
	for (unsigned int iV=0; iV<allDWIs.dwiVals.size(); iV++ ) dwiVals[iV] = log( (MyRealType)std::max( allDWIs.dwiVals[iV], m_LogMinArgumentValue ) );
	}

      if ( !m_DoWeightedLS )
	{
	sh::generateSHBasisMatrixPseudoInverse<MyRealType>( allDWIs.gradients, m_Lambda, m_NumTerms, m_UsedOrder, sh_basis_mat_pi );
	
	//generate intensity values
	
	s_values_new = m_sh_basis_mat_new * sh_basis_mat_pi * dwiVals;
	
	}
      else  // do the weighted least squares approximation with Huber
	// function instead
	{
	// need to make the matrices the correct size first
	
	unsigned int numoriggvectors = allDWIs.gradients.rows();
	
	B.set_size( numoriggvectors, m_NumTerms );
	sh::computeSHOrigBasisMat<MyRealType>( B, allDWIs.gradients, m_NumTerms );
	
	L.set_size( m_NumTerms );
	sh::computeSHRegularizationMatrix( L, m_NumTerms, m_UsedOrder );
	
	W.set_size( numoriggvectors );
	
	sh_basis_mat_pi.set_size( m_NumTerms, m_NumTerms );
	
	sh::generateSHBasisMatrixPseudoInversePrecomputed<MyRealType>( B, L, m_Lambda, sh_basis_mat_pi );
	// first compute a solution for the unweighted problem then do
	// some iterations for the weighting step
	
	vnl_vector<MyRealType> shCoeffs = sh_basis_mat_pi * dwiVals;
	
	for ( unsigned int iIter=0; iIter<m_NrOfWLSIterations; iIter++ )
	  {
	  sh::computeSHHuberWeightMatrix<MyRealType>( W, numoriggvectors, shCoeffs, B, dwiVals, m_HuberC, m_RiceSigma, nrOfDWIOutliers );
	  sh::generateWeightedSHBasisMatrixPseudoInversePrecomputed( B, L, W, m_Lambda, sh_basis_mat_pi );
	  shCoeffs = sh_basis_mat_pi*dwiVals;
	  }
	
	s_values_new = m_sh_basis_mat_new * sh_basis_mat_pi * dwiVals;
	
	}
      }
    else
      {
      // in case of a voxel which is not in the mask
      allBaselines.goodVoxel = false;
      allDWIs.goodVoxel = false;
      }

    // first write out the baselines

    itk::VariableLengthVector<DWIPixelType> outputVals(  m_numnewgvectors + m_NrOfBaselines );

    // first output the baselines

    // TODO: Need to have some form of average for the baselines

    if ( allBaselines.goodVoxel && isInMask )
      {
      // average them over all the cases using the geometric mean
      
      std::vector<MyRealType> averagedBaselines;
      computeAveragedBaselines( averagedBaselines, allBaselines, iNrOfBaselines, m_InterpolationType, m_AveragingType, nrOfBaselineOutliers );

      for ( unsigned int iV=0; iV<m_NrOfBaselines; iV++ )
	{
	outputVals[iV] = (DWIPixelType)round( m_ScalingFactor*averagedBaselines[iV] );
	}
      }
    else
      {
      for ( unsigned int iV=0; iV<m_NrOfBaselines; iV++ )
	{
	outputVals[iV] = (DWIPixelType)0;
	}
      }

    if ( !m_DoWeightedLS )
      {
      outlierit.Set( 0 );
      }
    else
      {
      if ( nrOfDWIOutliers + nrOfBaselineOutliers > 1000 )
	{
	std::cout << "Something is wrong here with the outlier detection = " << nrOfDWIOutliers + nrOfBaselineOutliers << std::endl;
	}
      outlierit.Set( nrOfDWIOutliers + nrOfBaselineOutliers );
      }

    if ( allDWIs.goodVoxel && isInMask )
      {
      if ( m_NoLogFit )
	{
	for ( unsigned int iV=0; iV<m_numnewgvectors; iV++ )
	  {
	  outputVals[iV+m_NrOfBaselines] = (DWIPixelType)round( m_ScalingFactor*s_values_new[iV] );
	  }
	}
      else
	{
	for ( unsigned int iV=0; iV<m_numnewgvectors; iV++ )
	  {
	  outputVals[iV+m_NrOfBaselines] = (DWIPixelType)round( m_ScalingFactor*exp( s_values_new[iV] ) );
	  }
	}
      }
    else
      {
      for ( unsigned int iV=0; iV<m_numnewgvectors; iV++ )
	{
	outputVals[iV+m_NrOfBaselines] = (DWIPixelType)0;
	}
      }

    // TODO: remove, for now just write the input to the output

    recondwiit.Set( outputVals );

    // increment all the iterators

    ++recondwiit;
    if ( recondwiit.IsAtEnd() ) noMoreElements = true;

    ++outlierit;
    if ( outlierit.IsAtEnd() ) noMoreElements = true;

    ++maskit;
    if ( maskit.IsAtEnd() ) noMoreElements = true;

    for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
      {
      ++(dwiits[iI]);
      if ( dwiits[iI].IsAtEnd() ) noMoreElements = true;
      }

    if ( !m_JustDoResampling )
      {
      for ( unsigned int iI=0; iI<nrOfDatasets; iI++ )
	{
	++(hfieldits[iI]);
	if ( hfieldits[iI].IsAtEnd() ) noMoreElements = true;
	}
      }

    if ( noMoreElements )
      {
      std::cout << "no more elements; thread id = " << threadId << std::endl;
      std::cout.flush();
      }

    progress.CompletedPixel();
    }

} // end update function  

template <class MyRealType, class DWIPixelType>
void 
DWIAtlasBuilder<MyRealType, DWIPixelType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Size : " << m_Size << std::endl;
  os << indent << "Origin: " << m_Origin << std::endl;
  os << indent << "Spacing: " << m_Spacing << std::endl;
}


} // end namespace itk

#endif
