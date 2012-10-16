/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToDTITubeSpatialObjectFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToDTITubeSpatialObjectFilter_h
#define __itkImageToDTITubeSpatialObjectFilter_h

namespace itk
{

/** \class ImageToDTITubeSpatialObjectFilter
 * \brief 
 *
 * ImageToDTITubeSpatialObjectFilter is the base class for all process objects that output
 * DTITubeSpatialObject data and require image data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ImageFilters
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
class ITK_EXPORT ImageToDTITubeSpatialObjectFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ImageToDTITubeSpatialObjectFilter  Self;
  typedef ProcessObject  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToDTITubeSpatialObjectFilter, DTITubeSpatialObjectSource);

  /** Create a valid output. */
  DataObject::Pointer  MakeOutput(unsigned int idx);

  /** Some Image related typedefs. */
  typedef   TInputImage                             InputImageType;
  typedef   typename InputImageType::Pointer        InputImagePointer;
  typedef   typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef   typename InputImageType::RegionType     InputImageRegionType; 
  typedef   typename InputImageType::PixelType      InputImagePixelType; 

  /** Some DTITubeSpatialObject related typedefs. */
  typedef   TOutputDTITubeSpatialObject                             OutputDTITubeSpatialObjectType;
  typedef   typename OutputDTITubeSpatialObjectType::Pointer        OutputDTITubeSpatialObjectPointer;

  /** Set the input image of this process object.  */
  void SetInput(unsigned int idx, const InputImageType *input);

  /** Get the input image of this process object.  */
  const InputImageType * GetInput(unsigned int idx) const;
  InputImageType * GetInput(unsigned int idx);

  /** Get the output DTITubeSpatialObject of this process object.  */
  const OutputDTITubeSpatialObjectType * GetOutput(void) const;
  OutputDTITubeSpatialObjectType * GetOutput(void);

  /** Prepare the output */
  void GenerateOutputInformation(void);
     
protected:
  ImageToDTITubeSpatialObjectFilter();
  ~ImageToDTITubeSpatialObjectFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
 
private:
  ImageToDTITubeSpatialObjectFilter(const ImageToDTITubeSpatialObjectFilter&); //purposely not implemented
  void operator=(const ImageToDTITubeSpatialObjectFilter&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToDTITubeSpatialObjectFilter.txx"
#endif

#endif
