/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorFileWriter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorFileWriter_h
#define __itkTensorFileWriter_h

#include <itkImageToImageFilter.h>
#include <itkDiffusionTensor3D.h>
#include <iostream>
#include <fstream>

#include <string>
#include <itkByteSwapper.h>

namespace itk
{

/** \class TensorFileWriter
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
template<typename T>
class ITK_EXPORT TensorFileWriter :
    public
ImageToImageFilter<Image<DiffusionTensor3D<T>, 3>,
                   Image<Vector<T,6>, 3> >
{
public:
//  typedef Image<Vector<typename TInputImage::PixelType::RealValueType,6>,3 > TOutputImage;
  typedef Image<DiffusionTensor3D<T>, 3> InputImageType;
  typedef Image<Vector<T,6>, 3> OutputImageType;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename InputPixelType::EigenVectorsMatrixType EigenVectorType;
  typedef typename InputPixelType::EigenValuesArrayType EigenValueType;
  typedef T RealType;


  /** Standard class typedefs. */
  typedef TensorFileWriter  Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType >  Superclass;

  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkGetMacro(FileName,std::string);
  itkSetMacro(FileName,std::string);
  
  /** Print internal ivars */
  void PrintSelf(std::ostream& os, Indent indent) const
    { this->Superclass::PrintSelf( os, indent ); }
  
  virtual void GenerateData();

protected:
  TensorFileWriter() {};
  virtual ~TensorFileWriter() {};

private:
  TensorFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string m_FileName;
};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorFileWriter.txx"
#endif


#endif
