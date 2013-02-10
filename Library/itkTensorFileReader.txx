#ifndef _itkTensorFileReader_txx
#define _itkTensorFileReader_txx
#include "itkTensorFileReader.h"

#include <itkByteSwapper.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <iostream>
#include <fstream>
#include <sstream>

namespace itk
{

template <class TOutputImage>
TensorFileReader<TOutputImage>
::TensorFileReader()
{
  this->m_FileName = "";
}

template <class TOutputImage>
TensorFileReader<TOutputImage>
::~TensorFileReader()
{

}

template <class TOutputImage>
void TensorFileReader<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  os << indent << "FileName: " << m_FileName << std::endl;
  os << indent << "Type:     " << m_Type << std::endl;
  os << indent << "NElem:    " << m_NElem << std::endl;
}

template <class TOutputImage>
void TensorFileReader<TOutputImage>
::GenerateOutputInformation()
{
  typename TOutputImage::Pointer output = this->GetOutput();

  if( m_FileName == "" )
    {
    throw ImageFileReaderException(__FILE__, __LINE__, "FileName must be specified");
    }

  SizeType      dimSize;
  double        spacing[TOutputImage::ImageDimension];
  double        origin[TOutputImage::ImageDimension];
  unsigned long nelem;

  char          line[256];
  std::ifstream instr(m_FileName.c_str() );
  std::string   buf;
  std::string   type;

  // Skip first four  lines, we assume we know what they mean
  instr.getline(line, 256);
  instr.getline(line, 256);
  instr.getline(line, 256);
  instr.getline(line, 256);

  // DIMESIONS size_x size_y size_z
  instr >> buf;
  instr >> dimSize[0];
  instr >> dimSize[1];
  instr >> dimSize[2];

  instr >> buf;
  instr >> spacing[0];
  instr >> spacing[1];
  instr >> spacing[2];

  instr >> buf;
  instr >> origin[0];
  instr >> origin[1];
  instr >> origin[2];

  instr >> buf;
  instr >> nelem;

  m_NElem = nelem;

  instr.getline(line, 256);
  instr.getline(line, 256);
  std::stringstream linestr(line);
  linestr >> buf;
  linestr >> buf;
  linestr >> type;

  m_Type = type;

  typename TOutputImage::DirectionType direction;
  for( int i = 0; i < 3; ++i )
    {
    for( int j = 0; j < 3; ++j )
      {
      direction[i][j] =  (i == j) ? 1 : 0;
      }
    }

  output->SetDirection( direction ); // Set the image direction cosines
  output->SetSpacing( spacing );     // Set the image spacing
  output->SetOrigin( origin );       // Set the image origin

  typedef typename TOutputImage::IndexType IndexType;

  IndexType start;
  start.Fill(0);

  ImageRegionType region;
  region.SetSize(dimSize);
  region.SetIndex(start);

  output->SetLargestPossibleRegion(region);

}

template <class TOutputImage>
void TensorFileReader<TOutputImage>
::GenerateData()
{
  std::ifstream instr(m_FileName.c_str() );
  // Skip the 9 header lines
  char line[256];

  for( int i = 0; i < 9; ++i )
    {
    instr.getline(line, 256);
    }

  typename TOutputImage::Pointer output = this->GetOutput();
  output->SetBufferedRegion( output->GetRequestedRegion() );

  output->Allocate();

  unsigned int elemsize = 8;
  int          type = 0;
  if( m_Type == "float" )
    {
    elemsize = 4;
    type = 0;
    }
  else if( m_Type == "double" )
    {
    elemsize = 8;
    type = 1;
    }
  else
    {
    std::cerr << "Unknown precision for tensor " << std::endl;
    }

  ImageRegionIteratorWithIndex<TOutputImage> it(output, output->GetRequestedRegion() );
  OutputPixelType                            pix;
  char*                                      buf = new char[9 * elemsize];
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    instr.read(buf, 9 * elemsize);
    if( type == 0 ) // Float
      {
      float* fldata = reinterpret_cast<float *>(buf);
      ByteSwapper<float>::SwapRangeFromSystemToBigEndian(fldata, 9);
      pix[0] = fldata[0];
      pix[1] = fldata[1];
      pix[2] = fldata[2];
      pix[3] = fldata[4];
      pix[4] = fldata[5];
      pix[5] = fldata[8];
      }
    else
      {
      double* ddata = reinterpret_cast<double *>(buf);
      ByteSwapper<double>::SwapRangeFromSystemToBigEndian(ddata, 9);
      pix[0] = ddata[0];
      pix[1] = ddata[1];
      pix[2] = ddata[2];
      pix[3] = ddata[4];
      pix[4] = ddata[5];
      pix[5] = ddata[8];
      }
    it.Set(pix);
    }
  delete[] buf;
}

}

#endif
