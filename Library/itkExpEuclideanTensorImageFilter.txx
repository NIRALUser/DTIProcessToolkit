// const double eps = 1e-16;
template <typename T>
void ExpEuclideanTensorImageFilter<T>::ThreadedGenerateData()
{
//  Superclass::GenerateInputRequestedRegion();
//  outputImage->SetRequrestedRegion(inputImage
  this->AllocateOutputs();

  const typename InputImageType::ConstPointer input(this->GetInput() );

  typename OutputImageType::Pointer output(this->GetOutput() );

  typename InputImageType::RegionType inputRequestedRegion(input->GetRequestedRegion() );

  typename OutputImageType::RegionType outputRequestedRegion(output->GetRequestedRegion() );

  ImageRegionConstIterator<InputImageType> it = ImageRegionConstIterator<InputImageType>(
      input, inputRequestedRegion);

  ImageRegionIterator<OutputImageType> oit = ImageRegionIterator<OutputImageType>(
      output, outputRequestedRegion);

  InputPixelType logvec;
  for( it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit )
    {
    logvec = it.Get();

    DiffusionTensor3D<double> tensor;
    for( int i = 0; i < 6; ++i )
      {
      tensor[i] = logvec[i];
      }
    tensor[1] /= sqrt(2.0);
    tensor[2] /= sqrt(2.0);
    tensor[4] /= sqrt(2.0);

    Vector<double, 3>    D;
    Matrix<double, 3, 3> U;

    tensor.ComputeEigenAnalysis(D, U);

    vnl_matrix_fixed<double, 3, 3> m;
    m.fill(0);
    m(0, 0) = exp(D[0]);
    m(1, 1) = exp(D[1]);
    m(2, 2) = exp(D[2]);

//    std::cout << U << std::endl;
//    std::cout << D << std::endl;

    vnl_matrix_fixed<double, 3, 3> res(U.GetVnlMatrix().transpose() * m * U.GetVnlMatrix() );

    OutputPixelType op;
    op[0] = res(0, 0);
    op[1] = res(0, 1);
    op[2] = res(0, 2);
    op[3] = res(1, 1);
    op[4] = res(1, 2);
    op[5] = res(2, 2);

    oit.Set(op);
    }

}
