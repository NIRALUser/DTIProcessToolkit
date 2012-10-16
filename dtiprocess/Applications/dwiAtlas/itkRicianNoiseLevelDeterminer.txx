/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRicianNoiseLevelDeterminer.txx,v $
  Language:  C++
  Date:      $Date: 2005-04-20 20:31:22 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkRicianNoiseLevelDeterminer_txx
#define _itkRicianNoiseLevelDeterminer_txx

#include "itkRicianNoiseLevelDeterminer.h"

#include "itkConstNeighborhoodIterator.h"
#include "vnl/vnl_math.h"


namespace itk {

template< class TImageType, class RealType >
RicianNoiseLevelDeterminer< TImageType, RealType >
::RicianNoiseLevelDeterminer() 
{
  m_RadiusEstimation.Fill(3);
  m_MinimumNumberOfUsedVoxelsEstimation = 1;
  m_MinimumNoiseSTD = 0;
  m_MaximumNoiseSTD = 200;
  m_HistogramResolutionFactor = 10;
  m_Verbose = false;
  m_HistogramFilename = "None";
}
    

 
template< class TImageType, class RealType >
void
RicianNoiseLevelDeterminer< TImageType, RealType >
::Compute( void )
{

  // first create a mask that contains all voxels different from zero

  typename ThresholdImageFilterType::Pointer zeroMaskImageFilter = ThresholdImageFilterType::New();

  zeroMaskImageFilter->SetInput( m_Input ); 
  zeroMaskImageFilter->ThresholdOutside(0,0);
  zeroMaskImageFilter->SetOutsideValue( 1 );

 // Compute the mean image for the input (which needs to be a baseline image)

  typedef itk::MaskedMeanImageFilter<ScalarImageType, ScalarRealImageType> MaskedMeanImageFilterType;

  typename MaskedMeanImageFilterType::Pointer maskedMeanImageFilter = MaskedMeanImageFilterType::New();
  maskedMeanImageFilter->SetInput( m_Input );
  maskedMeanImageFilter->SetMinimumNumberOfUsedVoxels( m_MinimumNumberOfUsedVoxelsEstimation );
  maskedMeanImageFilter->SetRadius( m_RadiusEstimation );

  // now generate the restricted mask, based on the mean of the filtered values and standard deviations
  
  typedef itk::LabelStatisticsImageFilter< ScalarRealImageType, ScalarImageType > LabelStatisticsImageFilterType;
  
  typename LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetInput( maskedMeanImageFilter->GetOutput() );
  labelStatisticsImageFilter->SetLabelInput( zeroMaskImageFilter->GetOutput() );
  labelStatisticsImageFilter->Update();

  RealType dMean = labelStatisticsImageFilter->GetMean( 1 );
  RealType dVariance = labelStatisticsImageFilter->GetVariance( 1 );
  RealType dSTD = sqrt(dVariance);
  RealType dMin = labelStatisticsImageFilter->GetMinimum( 1 );
  RealType dMax = labelStatisticsImageFilter->GetMaximum( 1 );

  if ( m_Verbose )
    {
    std::cout << "labelStatisticsImageFilter: " << std::endl;
    std::cout << "mean = " << dMean << std::endl;
    std::cout << "variance = " << dVariance << std::endl;
    std::cout << "standard deviation = " << dSTD << std::endl;
    std::cout << "min = " << dMin << std::endl;
    std::cout << "max = " << dMax << std::endl << std::endl;
    }
  // now set up the range we are interested in

  RealType dSigmaFac = 2.0;

  RealType dDesiredMin = dMean-dSigmaFac*dSTD-0.5;
  if ( dDesiredMin<dMin ) dDesiredMin = dMin-0.5;
  if ( dDesiredMin<=0 ) dDesiredMin = 0;

  RealType dDesiredMax = dMean+dSigmaFac*dSTD;
  if ( dDesiredMax>dMax ) dDesiredMax = dMax;
  
  // and extract values only within this range

  typedef itk::BinaryThresholdImageFilter<ScalarRealImageType,ScalarImageType> BinaryThresholdImageFilterType;

  typename BinaryThresholdImageFilterType::Pointer binaryRestrictedMask = BinaryThresholdImageFilterType::New();

  binaryRestrictedMask->SetInput( maskedMeanImageFilter->GetOutput() );
  binaryRestrictedMask->SetLowerThreshold( dDesiredMin );
  binaryRestrictedMask->SetUpperThreshold( dDesiredMax );
  binaryRestrictedMask->SetInsideValue( 1 );
  binaryRestrictedMask->SetOutsideValue( 0 );

  typedef itk::AndImageFilter<ScalarImageType,ScalarImageType,ScalarImageType> AndImageFilterType;

  typename AndImageFilterType::Pointer combinedMaskFilter = AndImageFilterType::New();

  combinedMaskFilter->SetInput1( binaryRestrictedMask->GetOutput() );
  combinedMaskFilter->SetInput2( zeroMaskImageFilter->GetOutput() );

  typedef itk::LabelStatisticsImageFilter< ScalarRealImageType, ScalarImageType > LabelStatisticsImageFilterType;

  typename LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilterRestricted = LabelStatisticsImageFilterType::New();

  labelStatisticsImageFilterRestricted->SetInput( maskedMeanImageFilter->GetOutput() );
  labelStatisticsImageFilterRestricted->SetLabelInput( combinedMaskFilter->GetOutput() );
  dLowerBound = dDesiredMin;
  dUpperBound = dDesiredMax;

  iNumBins = (int)round( (dUpperBound-dLowerBound)*m_HistogramResolutionFactor );
  
  if ( m_Verbose )
    {
    std::cout << "number of bins = " << iNumBins << std::endl;
    std::cout << "Histogram settings: iNumBins = " << iNumBins << "  dLowerBound = " << dLowerBound << "  dUpperBound = " << dUpperBound << std::endl;
    }

  labelStatisticsImageFilterRestricted->UseHistogramsOn();
  labelStatisticsImageFilterRestricted->SetHistogramParameters( iNumBins, dLowerBound, dUpperBound);
  labelStatisticsImageFilterRestricted->Update();

  RealType dMeanR = labelStatisticsImageFilterRestricted->GetMean( 1 );
  RealType dVarianceR = labelStatisticsImageFilterRestricted->GetVariance( 1 );
  RealType dSTDR = sqrt(dVarianceR);
  RealType dMinR = labelStatisticsImageFilterRestricted->GetMinimum( 1 );
  RealType dMaxR = labelStatisticsImageFilterRestricted->GetMaximum( 1 );

  if ( m_Verbose )
    {
    std::cout << "Restricted values:" << std::endl;
    std::cout << "mean = " << dMeanR << std::endl;
    std::cout << "variance = " << dVarianceR << std::endl;
    std::cout << "standard deviation = " << dSTDR << std::endl;
    std::cout << "min = " << dMinR << std::endl;
    std::cout << "max = " << dMaxR << std::endl;
    }

  // get the histogram

  typedef typename LabelStatisticsImageFilterType::HistogramPointer HistogramPointer;
  typedef typename LabelStatisticsImageFilterType::HistogramType HistogramType;

  HistogramPointer hp =    labelStatisticsImageFilterRestricted->GetHistogram( 1 );

  // iterate through this thing

  unsigned int iSize = (hp->GetSize())[0]; // this is a one-dimensional histogram

  int iCurrentMaxIndex = 0;
  int iCurrentMaxFrequency = 0;

  std::ofstream outputStream;
  std::string histogramFileName = "histogram.dat";

  if ( m_HistogramFilename.compare("None")!=0 ) 
    {
    outputStream.open( histogramFileName.c_str() );
    }  

  for ( unsigned int iI=0; iI<iSize; iI++ ) 
    {

    RealType dCurrentBinValue = (hp->GetMeasurementVector( iI ))[0];
    
    if ( (int)(hp->GetFrequency( iI ))>iCurrentMaxFrequency && (dCurrentBinValue<=m_MaximumNoiseSTD) && (dCurrentBinValue>=m_MinimumNoiseSTD) ) 
      {
      iCurrentMaxIndex = iI;
      iCurrentMaxFrequency = (int)hp->GetFrequency( iI );
      }

    if ( m_HistogramFilename.compare("None")!=0 ) 
      {
      outputStream << hp->GetFrequency( iI ) << " " << (hp->GetMeasurementVector( iI ))[0] << std::endl;
      }

    }

  if ( m_HistogramFilename.compare("None")!=0 ) 
    {
    outputStream.close();
    }

  /*std::cout << "max frequency = " << iCurrentMaxFrequency << " at index = " << iCurrentMaxIndex << " which corresponds to = " << hp->GetMeasurementVector( iCurrentMaxIndex ) << std::endl;

  std::cout << "dLowerBound = " << dLowerBound << std::endl;
  std::cout << "dUpperBound = " << dUpperBound << std::endl;
  std::cout << "dStep = " << dStep << std::endl;*/

  RealType dNoiseSTD_Estimated = sqrt(2/M_PI)*((hp->GetMeasurementVector( iCurrentMaxIndex ))[0]);

  if ( m_Verbose )
    {
    std::cout << "Estimated noise standard deviation is = " << dNoiseSTD_Estimated << std::endl; 
    }

  
  m_Output = dNoiseSTD_Estimated;

}
    
template< class TImageType, class RealType >
void
RicianNoiseLevelDeterminer< TImageType, RealType >
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius estimation: " << m_RadiusEstimation << std::endl;
  os << indent << "Histogram resolution factor: " << m_HistogramResolutionFactor << std::endl;
  os << indent << "Minimum number of used voxels estimation: " << m_MinimumNumberOfUsedVoxelsEstimation << std::endl;
  os << indent << "Histogram filename: " << m_HistogramFilename << std::endl;

}

    
} // end of namespace itk 


#endif
