/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDWIAtlasBuilder.h,v $
  Language:  C++
  Date:      $Date: 2007-01-25 09:18:58 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDWIAtlasBuilder_h
#define __itkDWIAtlasBuilder_h

#include "itkImageSource.h"
#include "itkConceptChecking.h"
#include "itkVectorImage.h"
#include "itkSimpleDataObjectDecorator.h"
#include "itkImageFileReader.h"
#include "itkCommand.h"

#include "auxFunctions.h"
#include "sphericalHarmonicsFunctions.h"

#include <vector>
#include <string>
#include <algorithm>

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_inverse.h>


#include "itkDeformationFieldJacobianFilter.h"
#include "itkDiffusionTensor3DReconstructionLinearImageFilter.h"
#include "itkRicianNoiseLevelDeterminer.h"
#include "itkExtractVolumeFilter.h"

namespace itk
{

// Currently not used
class ConsoleProgressCommand : public itk::Command 
{
public:
  typedef  ConsoleProgressCommand   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  ConsoleProgressCommand(): m_MaxProgress(100), m_NumTicks(50), m_Progress(0) { m_LastPercent = 0; firstRun = true; };
  
public:
  
  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }
  
  void SetMaxProgress( unsigned long mP ) 
  {
    m_MaxProgress = mP;
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    if( itk::StartEvent().CheckEvent( &event ) )
      {
      std::cout << "Running" << std::endl;
      std::cout << "0%                    50%                    100%" << std::endl;
      std::cout << "=================================================" << std::endl;
      

      }
    else if( itk::EndEvent().CheckEvent( &event )) 
      {
      std::cout << std::endl;
      }
    else if( itk::ProgressEvent().CheckEvent( &event )) 
      {
      ++m_Progress;
      unsigned int currentPerc = (unsigned int)round(100*((double)m_Progress)/m_MaxProgress);
      if ( (currentPerc%2==0 && (currentPerc>m_LastPercent)) | firstRun )
	{
        std::cout << "*";
	std::cout.flush();
	firstRun = false;
	}
      }
  }
private:
  unsigned long m_MaxProgress;
  unsigned long m_NumTicks;
  unsigned long m_Progress;
  unsigned int m_LastPercent;
  bool firstRun;
};


/** \class DWIAtlasBuilder
 *  \brief This class takes creates a DWI atlas form a set of DWI images.
 *
 */

const unsigned int DIM = 3;
const char* NRRD_MEASUREMENT_KEY = "NRRD_measurement frame";

template <class MyRealType, class DWIPixelType>
class ITK_EXPORT DWIAtlasBuilder :
  public ImageSource<VectorImage<  DWIPixelType , DIM > >
{
public:

  /** Standard class typedefs. */
  typedef DWIAtlasBuilder                Self;
  typedef ImageSource< VectorImage< DWIPixelType,  DIM > > Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>             ConstPointer;
  
  typedef VectorImage< DWIPixelType, DIM > OutputImageType; 
  typedef typename Superclass::Pointer    OutputImagePointer;
  typedef typename OutputImageType::SpacingType SpacingType;
  typedef typename OutputImageType::PointType   PointType;

  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(DWIAtlasBuilder,ImageSource);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Some convenient typedefs. */
  typedef VectorImage<DWIPixelType, DIM> VectorImageType;
  typedef Image<DWIPixelType,DIM> ScalarImageType;

  typedef ImageFileReader<VectorImageType> FileReaderType;
  typedef ImageFileReader<ScalarImageType> ScalarFileReaderType;

  typedef RicianNoiseLevelDeterminer< ScalarImageType, RealType > RicianNoiseLevelDeterminerType;
  typedef ExtractVolumeFilter< VectorImageType, ScalarImageType > ExtractInputVolumeFilterType;

  typedef DeformationFieldJacobianFilter<DeformationImageType, MyRealType> MyJacobianFilterType;

  typedef typename MyJacobianFilterType::OutputImageType MyJacobianImageType;
  typedef typename MyJacobianImageType::PixelType  MyJacobianType;

  typedef DiffusionTensor3DReconstructionLinearImageFilter<DWIPixelType, MyRealType>
DiffusionEstimationFilterType;
  typedef typename DiffusionEstimationFilterType::GradientDirectionContainerType GradientDirectionContainerType;

  typedef vnl_vector_fixed<MyRealType, DIM> MyGradientType;

  typedef typename VectorImageType::SizeType SizeType;

 /**  Command for observing progress of internal pipeline filters */
  typedef typename ConsoleProgressCommand::Pointer ProgressCommandPointer;

  /** Since a string is not a dataobject, we use the decorator to push
   *  it down the pipeline */ 
  typedef SimpleDataObjectDecorator< std::string*  > InputStringObjectType; 
  
  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, DIM );
  
  /** Set/Get the input of this process object.  */
  /* We have text files as input here for data that needs to be loaded */

  virtual void SetInput( const std::string* sCaseFile ); // TODO: Change

  /** Set the spacing (size of a pixel) of the image. 
   *  \sa GetSpacing() */
  itkSetMacro(Spacing,SpacingType);
  virtual void SetSpacing(const double* values);

  /** Get the spacing (size of a pixel) of the image. 
   * For ImageBase and Image, the default data spacing is unity. */
  itkGetConstReferenceMacro(Spacing,SpacingType);

  /** Set the origin of the image. 
   * \sa GetOrigin() */
  itkSetMacro(Origin,PointType);
  virtual void SetOrigin(const double* values);

  const InputStringObjectType* GetInput( void );

  ScalarImageType* GetOutlierImage() { return OutlierImage; };

  /** Set verbose mode. */
  itkSetMacro( Verbose, bool);
  itkGetMacro( Verbose, bool);
  itkBooleanMacro( Verbose );

  itkSetMacro( IsHField, bool );
  itkGetMacro( IsHField, bool );
  itkBooleanMacro( IsHField );

  itkSetMacro( JustDoResampling, bool );
  itkGetMacro( JustDoResampling, bool );
  itkBooleanMacro( JustDoResampling );

  itkSetMacro( NoLogFit, bool );
  itkGetMacro( NoLogFit, bool );
  itkBooleanMacro( NoLogFit );

  itkSetMacro( GradientVectorFile, std::string );
  itkGetMacro( GradientVectorFile, std::string );

  itkSetMacro( MaskImageFileName, std::string );
  itkGetMacro( MaskImageFileName, std::string );

  itkSetMacro( SHOrder, unsigned int );
  itkGetMacro( SHOrder, unsigned int );

  itkSetMacro( Lambda, MyRealType );
  itkGetMacro( Lambda, MyRealType );

  itkSetMacro( ScalingFactor, MyRealType );
  itkGetMacro( ScalingFactor, MyRealType );

  itkSetMacro( LogMinArgumentValue, DWIPixelType );
  itkGetMacro( LogMinArgumentValue, DWIPixelType );

  itkSetMacro( HuberC, MyRealType );
  itkGetMacro( HuberC, MyRealType );

  itkSetMacro( RiceSigma, MyRealType );
  itkGetMacro( RiceSigma, MyRealType );

  itkSetMacro( DoWeightedLS, bool );
  itkGetMacro( DoWeightedLS, bool );
  itkBooleanMacro( DoWeightedLS );

  itkSetMacro( NrOfWLSIterations, unsigned int );
  itkGetMacro( NrOfWLSIterations, unsigned int );

  void SetInterpolationType( std::string interpolationType );
  std::string GetInterpolationType();

  void SetAveragingType( std::string averagingType );
  std::string GetAveragingType();

  void SetNrOfThreads( int iNrOfThreads );
  unsigned int GetNrOfThreads();

 /** Get the origin of the image.  */
  itkGetConstReferenceMacro(Origin,PointType);

  
protected:
  DWIAtlasBuilder();
  ~DWIAtlasBuilder();

  typedef struct {
    std::vector<DWIPixelType> dwiVals;
    vnl_matrix<MyRealType> gradients;
    std::vector<MyRealType> interpolationWeights;
    std::vector<bool> isBaseline;
    unsigned int NDWI;
    unsigned int iNrOfBaselinesPerVolume;
    bool goodVoxel;
  } STransformedGradientInformationType;

  void LoadDataAndInitialize();

  vnl_matrix<MyRealType>
  rotateGradients( typename GradientDirectionContainerType::Pointer gradientContainer, itk::Matrix<MyRealType, 3, 3> jacobian, std::vector<bool>& isBaseline, unsigned int &iNrOfBaselines, bool bypassRotation = false );

  void getGradients( typename DiffusionEstimationFilterType::GradientDirectionContainerType & gradientContainer, itk::MetaDataDictionary & dict,  vnl_matrix<MyRealType> imgf, unsigned int& iNrOfBaselines, std::vector<unsigned int> &vecBaselineIndices, bool bVERBOSE=false );

  void extractDesiredBaselinesAndDWIs( STransformedGradientInformationType& allBaselines, STransformedGradientInformationType& allDWIs, const STransformedGradientInformationType* transformedInformation, unsigned int nrOfDataSets , unsigned int interpolationType, unsigned int averagingType );

  void computeAveragedBaselines( std::vector<MyRealType>& averagedBaselines, STransformedGradientInformationType& allBaselines, unsigned int iNrOfBaselines, unsigned int interpolationType, unsigned int averagingType, unsigned int & nrOfBaselineOutliers );

  virtual void GenerateOutputInformation();

  // threaded version to generate data
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId );
  void AfterThreadedGenerateData( void );
  void BeforeThreadedGenerateData( void );

  SizeType        m_Size;
  SpacingType     m_Spacing;
  PointType       m_Origin;
  
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:

  std::vector<std::string> dwiFiles;
  std::vector<std::string> deformationFiles;
  unsigned int nrOfDatasets;
  
  typename FileReaderType::Pointer *dwireader;
  DeformationImageType::Pointer *deformation;
  typename DiffusionEstimationFilterType::GradientDirectionContainerType::Pointer *gradientContainers;
  typename MyJacobianFilterType::Pointer *jacobian;

  typename ScalarFileReaderType::Pointer maskReader;

  typename ScalarImageType::Pointer OutlierImage;
  typename ScalarImageType::Pointer MaskImage;

  DWIAtlasBuilder(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_Verbose;
  bool m_IsHField;
  bool m_NoLogFit;
  bool m_DoWeightedLS;

  unsigned int m_SHOrder;
  std::string m_GradientVectorFile;
  std::string m_MaskImageFileName;

  bool m_JustDoResampling;
  MyRealType m_Lambda;

  MyRealType m_ScalingFactor; // factor to scale up reconstructed
			      // measurements before they are put into
			      // the original datatype again; this is
			      // useful to make full use of the range
			      // of the output data for mis-scaled
			      // data acquisitions

  DWIPixelType m_LogMinArgumentValue; // minimum value for logarithm
			      // computations; ideally there should
			      // never be any zeros in the
			      // computations, if so replace the
			      // values (only for the log computations)
			      // by this value (default 1) to make the
			      // transformation work

  int m_NumberOfThreads;

  unsigned int m_InterpolationType;
  unsigned int m_AveragingType;

  unsigned int m_UsedOrder;
  unsigned int m_NumTerms;

  unsigned int m_NrOfBaselines;  // the number of baselines
  unsigned int m_numnewgvectors; // the number of new gradients

  vnl_matrix<MyRealType> m_DesiredGradients;
  vnl_matrix<MyRealType> m_sh_basis_mat_new;


  // robust estimation parameters

  MyRealType m_HuberC;   // coefficient which specifies the switch
		       // between quadratic and magnitude loss
		       // function

  MyRealType m_RiceSigma;
  unsigned int m_NrOfWLSIterations;
  
  // end robust estimation parameters

  ProgressCommandPointer m_ConsoleProgressCommandPointer;

  const static unsigned int NCONTROLPOINTS = 8;

  enum {NEAREST_NEIGHBOR, LINEAR, USE_ALL_WITH_WEIGHTING};
  enum {ALGEBRAIC, GEOMETRIC};
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDWIAtlasBuilder.txx"
#endif

#endif
