#ifndef _itkTensorFileWriter_txx
#define _itkTensorFileWriter_txx
#include "itkTensorFileWriter.h"
#include <itkImageRegionConstIterator.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ios>

namespace itk{

template<typename T>
void
TensorFileWriter<T>::GenerateData()
{
  std::ofstream outstr(m_FileName.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
  typename InputImageType::ConstPointer in = this->GetInput(0);
  typename InputImageType::RegionType inregion = in->GetLargestPossibleRegion();
  outstr << std::setprecision(17) << std::scientific;
  outstr << "# vtk DataFile Version 2.0" << std::endl;
  outstr << "Diffusion Tensor Data" << std::endl;
  outstr << "BINARY" << std::endl;
  outstr << "DATASET STRUCTURED_POINTS" << std::endl;
  outstr << "DIMENSIONS " << inregion.GetSize()[0] << " " <<
    inregion.GetSize()[1] << " " << inregion.GetSize()[2] << std::endl;
  outstr << "SPACING " << in->GetSpacing()[0] << " " <<
    in->GetSpacing()[1] << " " << in->GetSpacing()[2] <<std::endl;
  outstr << "ORIGIN " << in->GetOrigin()[0] << " " <<
    in->GetOrigin()[1] << " " << in->GetOrigin()[2] <<std::endl;
  outstr << "POINT_DATA " << inregion.GetSize()[0] * inregion.GetSize()[1]  * inregion.GetSize()[2]  << std::endl;
  outstr << "TENSORS tensor float" << std::endl;


  typedef ImageRegionConstIterator<InputImageType> IteratorType;
  IteratorType it(in,inregion);
  for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    InputPixelType t = it.Get();
    float writedata[9];
    writedata[0] = t[0];
    writedata[1] = t[1];
    writedata[2] = t[2];
    writedata[3] = t[1];
    writedata[4] = t[3];
    writedata[5] = t[4];
    writedata[6] = t[2];
    writedata[7] = t[4];
    writedata[8] = t[5];

    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(writedata, 9);
    outstr.write(reinterpret_cast<const char*>(writedata),9*sizeof(float));

    }
}
}
#endif
