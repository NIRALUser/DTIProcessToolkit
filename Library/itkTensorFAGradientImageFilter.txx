#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkNthElementImageAdaptor.h>

//const double eps = 1e-16;
template<typename T>
void TensorFAGradientImageFilter<T>::GenerateData()
{
//  Superclass::GenerateInputRequestedRegion();
//  outputImage->SetRequrestedRegion(inputImage
  this->AllocateOutputs();

  typedef NthElementImageAdaptor<InputImageType, T>                                ElementAdaptorType;
  typedef typename ElementAdaptorType::Pointer                                     ElementAdaptorPointer;

  typedef OutputImageType ComponentGradientImageType;

  typedef itk::GradientRecursiveGaussianImageFilter<ElementAdaptorType,
                                                    ComponentGradientImageType>                                       TensorComponentGradientType;
  typedef typename TensorComponentGradientType::Pointer                            TensorComponentGradientPointer;

  std::vector<ElementAdaptorPointer> componentadaptors(6);
  std::vector<TensorComponentGradientPointer> componentgradientfilters(6);

  const typename InputImageType::Pointer input(const_cast<InputImageType*>(this->GetInput()));

  for(int i = 0; i < 6; ++i)
    {
    componentadaptors[i] = ElementAdaptorType::New();
    componentadaptors[i]->SetImage(input);
    componentadaptors[i]->SelectNthElement(i);
    componentadaptors[i]->Update();

    componentgradientfilters[i] = TensorComponentGradientType::New();
    componentgradientfilters[i]->SetInput(componentadaptors[i]);
    componentgradientfilters[i]->SetSigma(m_Sigma);
    componentgradientfilters[i]->Update();

    }

  typename OutputImageType::Pointer output(this->GetOutput());

  typename InputImageType::RegionType inputRequestedRegion(input->GetRequestedRegion());

  typename OutputImageType::RegionType outputRequestedRegion(output->GetRequestedRegion());

  ImageRegionConstIterator<InputImageType> it = ImageRegionConstIterator<InputImageType>(
    input, inputRequestedRegion);

  typedef ImageRegionConstIterator<ComponentGradientImageType> ComponentGradientIterator;

  std::vector<ComponentGradientIterator> git(6);

  for(int i = 0; i < 6; ++i)
    {
    git[i] = ComponentGradientIterator(componentgradientfilters[i]->GetOutput(),inputRequestedRegion);

    git[i].GoToBegin();
    }

  ImageRegionIterator<OutputImageType> oit = ImageRegionIterator<OutputImageType>(
    output, outputRequestedRegion);

  for(it.GoToBegin(), oit.GoToBegin();!it.IsAtEnd(); ++it, ++oit)
    {

    InputPixelType D = it.Get();

    T & Dxx = D[0];
    T & Dxy = D[1];
    T & Dxz = D[2];
    T & Dyy = D[3];
    T & Dyz = D[4];
    T & Dzz = D[5];

    T J2 = Dxx*Dyy + Dxx*Dzz + Dyy*Dzz - Dxy*Dxy - Dxz*Dxz - Dyz*Dyz;
    T J4 = Dxx*Dxx + Dyy*Dyy + Dzz*Dzz + 2*Dxy*Dxy + 2*Dxz*Dxz + 2*Dyz*Dyz;

    OutputPixelType gradDxx = git[0].Get();
    OutputPixelType gradDxy = git[1].Get();
    OutputPixelType gradDxz = git[2].Get();
    OutputPixelType gradDyy = git[3].Get();
    OutputPixelType gradDyz = git[4].Get();
    OutputPixelType gradDzz = git[5].Get();

    OutputPixelType gradJ2 =
      (Dyy+Dzz)*gradDxx +
      (Dxx+Dzz)*gradDyy +
      (Dxx+Dyy)*gradDzz -
      2*Dxy*gradDxy        -
      2*Dxz*gradDxz        -
      2*Dyz*gradDyz;

    OutputPixelType gradJ4 =
      2*Dxx*gradDxx        +
      2*Dyy*gradDyy        +
      2*Dzz*gradDzz        +
      4*Dxy*gradDxy        +
      4*Dxz*gradDxz        +
      4*Dyz*gradDyz;

    OutputPixelType op = (J2*gradJ4-J4*gradJ2)/(2*J4*J4*sqrt(1-J2/J4));
    oit.Set(op);

    for(int i = 0; i < 6; ++i)
      {
      ++git[i];
      }

    }

}
