/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarImageToGreyLevelCooccurrenceMatrixGenerator.h,v $
  Language:  C++
  Date:      $Date: 2005-08-24 15:16:54 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRicianNoiseLevelDeterminer_h
#define __itkRicianNoiseLevelDeterminer_h

#include "itkImage.h"
#include "itkVectorContainer.h"
#include "itkObject.h"
#include "itkNumericTraits.h"
#include "itkMacro.h"

#include "itkMaskedMeanImageFilter.h"

#include <itkThresholdImageFilter.h>
#include <itkThresholdLabelerImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

#include <itkAndImageFilter.h>

#include <fstream>

namespace itk {

/** \class RicianNoiseLevelDeterminer 
*  \brief Computes the Rician noise level given a baseline image
*
* Author: Marc Niethammer
*/
    
template< class TImageType, class RealType >
class RicianNoiseLevelDeterminer : public Object
{
public:
  /** Standard typedefs */
  typedef RicianNoiseLevelDeterminer Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImageType::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RicianNoiseLevelDeterminer, Object);
  
  /** standard New() method support */
  itkNewMacro(Self) ;

  typedef TImageType ScalarImageType;
  typedef typename itk::Image<RealType,ImageDimension> ScalarRealImageType;  

  typedef typename ScalarImageType::ConstPointer                ImageConstPointer;
  typedef typename ScalarImageType::SizeType                    RadiusType;

  typedef itk::ThresholdImageFilter<ScalarImageType> ThresholdImageFilterType;
  
  /** Triggers the Computation of Rician noise level */
  void Compute( void );
  
  /** Connects the input image for which the histogram is going to be computed */
  itkSetConstObjectMacro( Input, ScalarImageType );
  itkGetConstObjectMacro( Input, ScalarImageType );

  // some macros

  void SetRadiusEstimation( unsigned int iRE ) { m_RadiusEstimation.Fill( iRE ); };
  itkGetMacro( RadiusEstimation, RadiusType );

  itkSetMacro( MinimumNoiseSTD, RealType );
  itkGetMacro( MinimumNoiseSTD, RealType );

  itkSetMacro( MaximumNoiseSTD, RealType );
  itkGetMacro( MaximumNoiseSTD, RealType );
  
  itkSetMacro( MinimumNumberOfUsedVoxelsEstimation, unsigned int );
  itkGetMacro( MinimumNumberOfUsedVoxelsEstimation, unsigned int );

  itkSetMacro( HistogramResolutionFactor, RealType );
  itkGetMacro( HistogramResolutionFactor, RealType );

  itkSetMacro( Verbose, bool );
  itkGetMacro( Verbose, bool );
  itkBooleanMacro( Verbose );

  itkSetMacro( HistogramFilename, std::string );
  itkGetMacro( HistogramFilename, std::string );

  /** Return the noise level
      \warning This output is only valid after the Compute() method has been invoked 
      \sa Compute */
  RealType GetOutput() { return m_Output; };
  
protected:
  RicianNoiseLevelDeterminer();
  virtual ~RicianNoiseLevelDeterminer() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
private:
  
  ImageConstPointer        m_Input;
  RealType                 m_Output;

  unsigned int m_MinimumNumberOfUsedVoxelsEstimation;

  RealType m_MinimumNoiseSTD;
  RealType m_MaximumNoiseSTD;

  RealType dLowerBound;
  RealType dUpperBound;
  int iNumBins;

  std::string m_HistogramFilename;

  RealType m_HistogramResolutionFactor;

  bool m_Verbose;

  RadiusType m_RadiusEstimation;
};


} // end of namespace itk 

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRicianNoiseLevelDeterminer.txx"
#endif

#endif
