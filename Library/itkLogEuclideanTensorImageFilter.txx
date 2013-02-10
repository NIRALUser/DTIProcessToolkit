// const double eps = 1e-16;
template <typename T>
void LogEuclideanTensorImageFilter<T>::ThreadedGenerateData()
{
//  Superclass::GenerateInputRequestedRegion();
//  outputImage->SetRequrestedRegion(inputImage
  this->AllocateOutputs();

  const typename InputImageType::ConstPointer input(this->GetInput() );

  typename OutputImageType::Pointer output(this->GetOutput() );

  typename InputImageType::RegionType inputRequestedRegion(input->GetRequestedRegion() );

  typename OutputImageType::RegionType outputRequestedRegion(output->GetRequestedRegion() );

  ImageRegionConstIteratorWithIndex<InputImageType> it(
    input, inputRequestedRegion);

  ImageRegionIterator<OutputImageType> oit = ImageRegionIterator<OutputImageType>(
      output, outputRequestedRegion);

  InputPixelType tensor;
  for( it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit )
    {
    tensor = it.Get();

    Vector<double, 3>    D;
    Matrix<double, 3, 3> U;

//    std::cout << it.GetIndex() << std::endl;
    DiffusionTensor3D<double> temp;
    for( int i = 0; i < 6; i++ )
      {
      temp[i] = tensor[i];
      }
    temp.ComputeEigenAnalysis(D, U);

    vnl_matrix_fixed<double, 3, 3> m;
    m.fill(0);
    m(0, 0) = D[0] > 0 ? log(D[0]) : -10;
    m(1, 1) = D[1] > 0 ? log(D[1]) : -10;
    m(2, 2) = D[2] > 0 ? log(D[2]) : -10;

//    std::cout << U << std::endl;
//    std::cout << D << std::endl;

    vnl_matrix_fixed<double, 3, 3> res(U.GetVnlMatrix().transpose() * m * U.GetVnlMatrix() );

    OutputPixelType op;
    op[0] = res(0, 0);
    op[1] = res(0, 1) * sqrt(2.0);
    op[2] = res(0, 2) * sqrt(2.0);
    op[3] = res(1, 1);
    op[4] = res(1, 2) * sqrt(2.0);
    op[5] = res(2, 2);

    oit.Set(op);
    }

}
