
#include <itkImageRegionConstIteratorWithIndex.h>

#include "itkImageToDTIStreamlineTractographyFilter.h"
#include "itkTensorPrincipalEigenvectorImageFilter.h"

namespace itk
{

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::ImageToDTIStreamlineTractographyFilter()
  : m_StepSize(0.5), m_MinimumFractionalAnisotropy(0.2), m_MaximumDirectionDeviation(1.0/1.414213),
    m_SourceLabel(2), m_TargetLabel(1)
{
  this->SetNumberOfRequiredInputs(2);

  m_TubeGroup = OutputGroupSpatialObjectType::New();
  m_TensorInterpolator = TensorInterpolateType::New();
  m_ROIInterpolator = ROIInterpolateType::New();
}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
void
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::SetTensorImage(const TTensorImage* timage)
{
  this->ProcessObject::SetNthInput(0, const_cast<TTensorImage *>( timage ));

}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
void
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::SetROIImage(const TROIImage* roiimage)
{
  this->ProcessObject::SetNthInput(1, const_cast<TROIImage *>( roiimage ));
  
}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
const TTensorImage*
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::GetTensorImage() const
{
  return this->GetInput(0);

}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
const TROIImage*
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::GetROIImage() const
{
  return dynamic_cast<const TROIImage*>(this->ProcessObject::GetInput(1));
}


template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
void
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::GenerateData()
{
  // Preprocessing: 
  this->PreprocessTensorImage();

  // Loop over ROI image
  typedef ImageRegionConstIteratorWithIndex<TensorImageType> TensorIteratorType;
  typedef ImageRegionConstIterator<ROIImageType> ROIIteratorType;

  TensorIteratorType tensorit(this->GetTensorImage(),this->GetTensorImage()->GetLargestPossibleRegion());
  ROIIteratorType roiit(this->GetROIImage(),this->GetROIImage()->GetLargestPossibleRegion());

  for(tensorit.GoToBegin(), roiit.GoToBegin();
      !tensorit.IsAtEnd();
      ++tensorit, ++roiit)
  {
    // For each pixel which is a source region
    if(roiit.Get() == m_SourceLabel)
    {
      IndexType ind = tensorit.GetIndex();
      itkDebugMacro(<<"Initalizing fiber from "<< ind);

      PointType pt;
      this->GetTensorImage()->TransformIndexToPhysicalPoint(ind,pt);

      // Find two initial starting directions
      TensorType tens = tensorit.Get();
      itkDebugMacro(<<"Pretensor: " << tens);

      EigenVectorType veca = Functor::TensorPrincipalEigenvectorFunction<TensorType, double>()(tens);
      EigenVectorType vecb = -veca;
      
      typename DTITubeSpatialObjectType::Pointer fiba, fibb, fib;

      // Track in first direction
      fiba = this->TrackFromPoint(pt,veca);

      // Track in second direction
      fibb = this->TrackFromPoint(pt,vecb);

      fib = DTITubeSpatialObjectType::New();
      typedef DTITubeSpatialObjectType::PointListType PointListType;
      PointListType newpoints(fiba->GetNumberOfPoints() + fibb->GetNumberOfPoints() - 1);
      std::copy(fiba->GetPoints().rbegin(), fiba->GetPoints().rend(),
                newpoints.begin());
      //  Plus one avoids double entering the start point
      std::copy(fibb->GetPoints().begin() + 1, fibb->GetPoints().end(),
                (newpoints.begin() + fiba->GetNumberOfPoints()));
      fib->SetPoints(newpoints);

      // Need to set spacing
      fib->SetSpacing(this->GetROIImage()->GetSpacing().GetDataPointer());

      // TODO: iterate over fiber and test if passed through target
      // region
      bool keep = false;
      if(!m_TargetLabel)
      {
        keep = true;
      }
      else
      {
        for(PointListType::const_iterator it = newpoints.begin();
            it != newpoints.end(); ++it)
        {
          ContinuousIndex<double, 3> cind;
          cind[0] = it->GetPosition()[0];
          cind[1] = it->GetPosition()[1];
          cind[2] = it->GetPosition()[2];
          PointType ptest;
          this->GetTensorImage()->TransformContinuousIndexToPhysicalPoint(cind, ptest);
          
          if(m_ROIInterpolator->Evaluate(ptest) == m_TargetLabel)
          {
            keep = true;
            break;
          }
        }
      }
      if(keep)
      {
        m_TubeGroup->AddSpatialObject(fib);
      }
    } // end if in source region

  } /// end loop over pixels
  m_TubeGroup->ComputeObjectToWorldTransform();
}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
typename DTITubeSpatialObject<3>::Pointer
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::TrackFromPoint(PointType pt,
                 EigenVectorType vec) const
{
  itkDebugMacro(<<"Tracking from " << pt << " in direction " << vec);

  std::vector<DTITubeSpatialObjectPointType> pointlist;
  DTITubeSpatialObjectTypePointer tube = DTITubeSpatialObjectType::New();
  
  bool stoppingcond = false;

  EigenVectorType nextvec;
  
  DTITubeSpatialObjectPointType tubept;
  ContinuousIndex<double, 3> nextpointind;
  
  this->GetTensorImage()->TransformPhysicalPointToContinuousIndex(pt, nextpointind);
  // Since this point came from an ROI label it ought to be guaranteed
  // to be in the image buffer.
  assert(m_TensorInterpolator->IsInsideBuffer(pt));
  tubept.SetPosition(nextpointind[0], nextpointind[1], nextpointind[2]);
  TensorType t = m_TensorInterpolator->Evaluate(pt);
  tubept.SetRadius(0.5);
  tubept.SetTensorMatrix(t);
  pointlist.push_back(tubept);

  do
  {
    PointType nextpt;
    try
    {
      nextpt = this->IntegrateOneStep(pt, vec, m_StepSize);
      nextvec = nextpt - pt;
    }
    catch (itk::ExceptionObject e)
    {
      break;
    }

    itkDebugMacro(<<"Next point: " << nextpt);

    typedef typename TensorInterpolateType::OutputType ArrayType;
    TensorType t = m_TensorInterpolator->Evaluate(nextpt);
    
    if(t.GetFractionalAnisotropy() < m_MinimumFractionalAnisotropy || 
       dot_product(nextvec.GetVnlVector().normalize(), vec.GetVnlVector().normalize()) < m_MaximumDirectionDeviation) // extremely arbitrarty
    {
      stoppingcond = true;
    }
    else
    {
      this->GetTensorImage()->TransformPhysicalPointToContinuousIndex(nextpt, nextpointind);

      tubept.SetPosition(nextpointind[0], nextpointind[1], nextpointind[2]);
      tubept.SetTensorMatrix(t);
      pointlist.push_back(tubept);
      pt = nextpt;
      vec = nextvec;
    }
    
  }
  while(!stoppingcond);

  tube->SetPoints(pointlist);

  return tube;
}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
typename TTensorImage::PointType
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::IntegrateOneStep(const PointType& pt,
                   const EigenVectorType& vec,
                   double h) const
{
  // Implement rk4
  EigenVectorType k1,k2,k3,k4;
  //EigenVectorType k1;

  PointType testpoint;

  // Evaluate next point using 4-order runge-kutta integrationn
  k1 = this->EvaluatePrincipalDiffusionDirectionAt(pt, vec);

  for(unsigned int i = 0; i < PointType::Dimension; ++i)
  {
    testpoint[i] = pt[i] + m_StepSize*k1[i]/2;
  }
  if(!m_TensorInterpolator->IsInsideBuffer(testpoint))
    throw itk::ExceptionObject("Tracked outside image buffer");

  k2 = this->EvaluatePrincipalDiffusionDirectionAt(testpoint, vec);

  for(unsigned int i = 0; i < PointType::Dimension; ++i)
  {
    testpoint[i] = pt[i] + m_StepSize*k2[i]/2;
  }
  if(!m_TensorInterpolator->IsInsideBuffer(testpoint))
    throw itk::ExceptionObject("Tracked outside image buffer");

  k3 = this->EvaluatePrincipalDiffusionDirectionAt(testpoint, vec);

  for(unsigned int i = 0; i < PointType::Dimension; ++i)
  {
    testpoint[i] = pt[i] + m_StepSize*k3[i];
  }
  if(!m_TensorInterpolator->IsInsideBuffer(testpoint))
    throw itk::ExceptionObject("Tracked outside image buffer");

  k4 = this->EvaluatePrincipalDiffusionDirectionAt(testpoint, vec);
  
  PointType geompoint;
  for(unsigned int i = 0; i < PointType::Dimension; ++i)
  {
    geompoint[i] = pt[i] + k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6;
//    geompoint[i] = pt[i] + h*k1[i];
  }

  if(!m_TensorInterpolator->IsInsideBuffer(geompoint))
    throw itk::ExceptionObject("Tracked outside image buffer");

  return geompoint;
}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
CovariantVector<double,3>
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::EvaluatePrincipalDiffusionDirectionAt(const PointType& pt,
                                        const EigenVectorType& vec) const
{
  TensorType tens = m_TensorInterpolator->Evaluate(pt);
  
  EigenVectorType pdd = Functor::TensorPrincipalEigenvectorFunction<TensorType, double>()(tens);

  itkDebugMacro(<<"PDD: " << pdd);

  if(pdd[0]*vec[0] + pdd[1]*vec[1] + pdd[2]*vec[2] < 0)
  {
    pdd = -pdd;
  }

  return pdd;
}

template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
void
ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
::PreprocessTensorImage() // BeforeThreadedGenerateData
{
  m_TensorInterpolator->SetInputImage(this->GetTensorImage());
  m_ROIInterpolator->SetInputImage(this->GetROIImage());

  // TODO: this should not be here
  m_TubeGroup->SetSpacing(this->GetTensorImage()->GetSpacing().GetDataPointer());
  m_TubeGroup->GetObjectToParentTransform()->SetOffset(this->GetTensorImage()->GetOrigin().GetDataPointer());
  
}

// template <class TTensorImage, class TROIImage, class TOutputSpatialObject>
// void ImageToDTIStreamlineTractographyFilter<TTensorImage,TROIImage,TOutputSpatialObject>
// ::GenerateOutputInformation()
// {
//   std::cout << "Generate Output Info" << std::endl;
// }

} // end namespace itk
