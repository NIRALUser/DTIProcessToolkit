#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

namespace itk
{

//const double eps = 1e-16;
template <typename TInputImage, typename TOutputImage>
void HFieldToDeformationFieldImageFilter<TInputImage,TOutputImage>::GenerateData()
{
//  Superclass::GenerateInputRequestedRegion();
//  outputImage->SetRequrestedRegion(inputImage
  this->AllocateOutputs();

  const typename InputImageType::ConstPointer input(this->GetInput() );

  typename OutputImageType::Pointer output(this->GetOutput());
  
  typename InputImageType::RegionType inputRequestedRegion(input->GetRequestedRegion());

  typename OutputImageType::RegionType outputRequestedRegion(output->GetRequestedRegion());

  ImageRegionConstIteratorWithIndex<InputImageType> it = ImageRegionConstIteratorWithIndex<InputImageType>(
    input, inputRequestedRegion);
  
  ImageRegionIterator<OutputImageType> oit = ImageRegionIterator<OutputImageType>(
    output, outputRequestedRegion);
  
  //  typename InputImageType::SpacingType spacing = input->GetSpacing();
  for(it.GoToBegin(), oit.GoToBegin();!it.IsAtEnd(); ++it, ++oit)
    {
    InputPixelType hvec = it.Get();
    typename OutputImageType::IndexType index = it.GetIndex();

    oit.Set(this->ComputeDisplacement(input,index,hvec));
    }

}

template <typename TInputImage, typename TOutputImage>
typename TOutputImage::PixelType 
HFieldToDeformationFieldImageFilter<TInputImage,TOutputImage>::ComputeDisplacement(typename InputImageType::ConstPointer input,
                                                                                   typename InputImageType::IndexType ind,
                                                                                   typename InputImageType::PixelType hvec)
{
  //typedef typename InputImageType::PixelType InputPixelType;
  typedef typename InputPixelType::ValueType CoordRepType;
  const unsigned int Dimension = InputImageType::ImageDimension;

  typedef itk::Point<CoordRepType, Dimension> PointType;
  PointType ipt, hpt;


  typedef ContinuousIndex<CoordRepType, Dimension> ContinuousIndexType;
  ContinuousIndexType hind;
  for(unsigned int i = 0; i < Dimension; ++i)
    {
    hind[i] = hvec[i];
    }

  input->TransformIndexToPhysicalPoint(ind,ipt);
  input->TransformContinuousIndexToPhysicalPoint(hind,hpt);

  return hpt - ipt;
}

} // namespace itk
