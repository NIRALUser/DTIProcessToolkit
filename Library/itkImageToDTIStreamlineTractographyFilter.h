#ifndef __itkImageToDTIStreamlineTractographyFilter_h
#define __itkImageToDTIStreamlineTractographyFilter_h

// Base class
#include "itkImageToDTITubeSpatialObjectFilter.h"

#include <itkGroupSpatialObject.h>
#include <itkDTITubeSpatialObject.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkTensorLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

namespace itk
{

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
class ITK_EXPORT ImageToDTIStreamlineTractographyFilter :
  public         ImageToDTITubeSpatialObjectFilter<TTensorImage, TOutputSpatialObject>
{
public:
  /** Standard class typedefs. */
  typedef ImageToDTIStreamlineTractographyFilter Self;
  typedef ProcessObject                          Superclass;
  typedef SmartPointer<Self>                     Pointer;
  typedef SmartPointer<const Self>               ConstPointer;

  typedef TOutputSpatialObject OutputGroupSpatialObjectType;
  typedef typename
    OutputGroupSpatialObjectType::Pointer       OutputGroupSpatialObjectPointer;

  // Image Types

  typedef TTensorImage                        TensorImageType;
  typedef typename TensorImageType::PixelType TensorType;
  typedef itk::CovariantVector<double, 3>     EigenVectorType;
  typedef TROIImage                           ROIImageType;
  typedef typename ROIImageType::PixelType    ROIPixelType;

  // Point and index types
  typedef typename TensorImageType::IndexType IndexType;
  typedef typename TensorImageType::PointType PointType;

  typedef TensorLinearInterpolateImageFunction<TensorImageType, double> TensorInterpolateType;
  typedef typename TensorInterpolateType::Pointer                       TensorInterpolatePointer;
  typedef NearestNeighborInterpolateImageFunction<ROIImageType, double> ROIInterpolateType;
  typedef typename ROIInterpolateType::Pointer                          ROIInterpolatePointer;

  // Spatial object Types
  typedef DTITubeSpatialObject<3>                   DTITubeSpatialObjectType;
  typedef typename DTITubeSpatialObject<3>::Pointer DTITubeSpatialObjectTypePointer;

  typedef DTITubeSpatialObjectPoint<3> DTITubeSpatialObjectPointType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  itkTypeMacro(ImageToDTIStreamlineTractographyFilter, ImageToDTITubeSpatialObjectFilter);

  itkGetMacro( StepSize, double );
  itkSetMacro( StepSize, double );

  itkGetMacro( MinimumFractionalAnisotropy, double );
  itkSetMacro( MinimumFractionalAnisotropy, double );

  itkGetMacro( MaximumAngleChange, double );
  itkSetMacro( MaximumAngleChange, double );

  itkGetMacro( SourceLabel, ROIPixelType );
  itkSetMacro( SourceLabel, ROIPixelType );

  itkGetMacro( TargetLabel, ROIPixelType );
  itkSetMacro( TargetLabel, ROIPixelType );

  itkGetMacro( ForbiddenLabel, ROIPixelType );
  itkSetMacro( ForbiddenLabel, ROIPixelType );

  itkGetMacro( WholeBrain, bool );
  itkSetMacro( WholeBrain, bool );
  itkBooleanMacro( WholeBrain );

  virtual void SetTensorImage(const TTensorImage* timage);

  virtual void SetROIImage(const TROIImage* roiimage);

  virtual const TTensorImage * GetTensorImage() const;

  virtual const TROIImage * GetROIImage() const;

  virtual void Update() ITK_OVERRIDE
  {
    this->GenerateData();
  }

  virtual OutputGroupSpatialObjectType * GetOutput()
  {
    return m_TubeGroup.GetPointer();
  }

protected:

  virtual void GenerateData() ITK_OVERRIDE;

  virtual void GenerateOutputInformation() ITK_OVERRIDE
  {
  };                                           // do nothing

  virtual DTITubeSpatialObjectTypePointer TrackFromPoint(PointType pt, EigenVectorType vec) const;

  virtual PointType IntegrateOneStep(const PointType& pt, const EigenVectorType& vec, double stepsize) const;

  // Preprocess tensor field to extract necessary information
  virtual void PreprocessTensorImage();

  virtual EigenVectorType EvaluatePrincipalDiffusionDirectionAt(const PointType& pt, const EigenVectorType& vec) const;

  ImageToDTIStreamlineTractographyFilter();
  virtual ~ImageToDTIStreamlineTractographyFilter()
  {
  };
private:
  ImageToDTIStreamlineTractographyFilter(const Self &); // purposely
  // not
  // implemented
  void operator=(const Self &); // purposely not implemented

  double m_StepSize;
  double m_MinimumFractionalAnisotropy;
  double m_MaximumAngleChange;

  ROIPixelType m_SourceLabel;
  ROIPixelType m_TargetLabel;
  ROIPixelType m_ForbiddenLabel;
  bool         m_WholeBrain;

  TensorInterpolatePointer m_TensorInterpolator;
  ROIInterpolatePointer    m_ROIInterpolator;

  OutputGroupSpatialObjectPointer m_TubeGroup;

}; // end class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToDTIStreamlineTractographyFilter.txx"
#endif

#endif
