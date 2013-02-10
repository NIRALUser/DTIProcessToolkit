// const double eps = 1e-16;
template <typename T>
void TensorFAHessianImageFilter<T>::GenerateData()
{
//  Superclass::GenerateInputRequestedRegion();
//  outputImage->SetRequrestedRegion(inputImage
  this->AllocateOutputs();

  typedef itk::NthElementImageAdaptor<InputImageType, T> ElementAdaptorType;
  typedef typename ElementAdaptorType::Pointer           ElementAdaptorPointer;

  typedef itk::GradientRecursiveGaussianImageFilter<TensorComponentAdaptorType,
                                                    ComponentGradientImageType>
    TensorComponentGradientType;
  typedef typename TensorComponentGradientType::Pointer TensorComponentGradientPointer;
  typedef itk::HessianRecursiveGaussianImageFilter<TensorComponentAdaptorType,
                                                   ComponentHessianImageType>
    TensorComponentHessianType;
  typedef typename TensorComponentHessianType::Pointer TensorComponentHessianPointer;

  std::vector<TensorComponentAdaptorPointer>  componentadaptors(6);
  std::vector<TensorComponentGradientPointer> componentgradientfilters(6);
  std::vector<TensorComponentHessianPointer>  componenthessianfilters(6);

  const typename InputImageType::ConstPointer input(this->GetInput() );
  for( int i = 0; i < 6; ++i )
    {
    componentadaptors[i] = NthElementImageAdaptor::New();
    componentadaptors[i]->SelectNthElement(i);
    componentadpators[i]->SetInput(input);

    componentgradientfilters[i]->SetInput(componentadaptors[i]->GetOutput() );
    componentgradientfilters[i]->SetSigma(m_Sigma);
    componentgradientfilters[i]->Update();

    componenthessianfilters[i]->SetInput(componentadaptors[i]->GetOutput() );
    componenthessianfilters[i]->SetSigma(m_Sigma);
    componenthessianfilters[i]->Update();

    }

  typename OutputImageType::Pointer output(this->GetOutput() );

  typename InputImageType::RegionType inputRequestedRegion(input->GetRequestedRegion() );

  typename OutputImageType::RegionType outputRequestedRegion(output->GetRequestedRegion() );

  ImageRegionConstIterator<InputImageType> it = ImageRegionConstIterator<InputImageType>(
      input, inputRequestedRegion);

  typedef ImageRegionConstIterator<ComponentGradientImageType> ComponentGradientIterator;

  typedef ImageRegionConstIterator<ComponentHessianImageType> ComponentHessianIterator;

  std::vector<ComponentGradientIterator> git(6);
  std::vector<ComponentHessianIterator>  hit(6);
  for( int i = 0; i < 6; ++i )
    {
    git = ComponentGradientIterator(componentgradientfilters[i]->GetOutput(), inputRequestedRegion);
    hit = ComponentHessianIterator(componenthessianfilters[i]->GetOutput(), inputRequestedRegion);

    git.GoToBegin();
    hit.GoToBegin();
    }

  ImageRegionIterator<OutputImageType> oit = ImageRegionIterator<OutputImageType>(
      output, outputRequestedRegion);
  for( it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit )
    {
    for( int i = 0; i < 6; ++i )
      {
      ++git[i];
      ++hit[i];
      }

    }

}
