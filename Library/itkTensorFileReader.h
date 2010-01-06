/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorFileReader.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorFileReader_h
#define __itkTensorFileReader_h

#include <itkImageSource.h>

#include <string>

namespace itk
{



/** \class TensorFileReader
 * \brief Computes the 6-element vector field that are the unique
 * elements of the matrix-logarithm of the tensor field.
 *
 * LogEuclideanImageFilter applies pixel-wise the invokation for
 * computing the matrix logarithm of every pixel. 
 * 
 * \sa DiffusionTensor3D
 * 
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 *
 */
template <class TOutputImage>
class ITK_EXPORT TensorFileReader :
    public ImageSource<TOutputImage>

{
public:
//  typedef Image<Vector<typename TInputImage::PixelType::RealValueType,6>,3 > TOutputImage;
  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename OutputPixelType::RealValueType RealType;

  typedef typename TOutputImage::SizeType  SizeType;

  /** The region of the output image. */
  typedef typename TOutputImage::RegionType  ImageRegionType;

  /** Standard class typedefs. */
  typedef TensorFileReader  Self;
  typedef ImageSource<OutputImageType >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkGetMacro(FileName,std::string);
  itkSetMacro(FileName,std::string);

  itkGetMacro(Type,std::string);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  virtual void GenerateOutputInformation(void);

protected:
  TensorFileReader();
  virtual ~TensorFileReader();

  virtual void GenerateData();

private:
  TensorFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string m_FileName;
  std::string m_Type;

  unsigned long m_NElem;

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorFileReader.txx"
#endif


#endif
