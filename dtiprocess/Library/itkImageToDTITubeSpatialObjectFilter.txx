/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToDTITubeSpatialObjectFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkImageToDTITubeSpatialObjectFilter_txx
#define _itkImageToDTITubeSpatialObjectFilter_txx
#include "itkImageToDTITubeSpatialObjectFilter.h"


namespace itk
{

/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::ImageToDTITubeSpatialObjectFilter()
{

  this->ProcessObject::SetNumberOfRequiredInputs(1);

  OutputDTITubeSpatialObjectPointer output
    = dynamic_cast<OutputDTITubeSpatialObjectType*>(this->MakeOutput(0).GetPointer()); 

  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());

}

/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::~ImageToDTITubeSpatialObjectFilter()
{
}
  

/**
 *   Make Ouput
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
DataObject::Pointer
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::MakeOutput(unsigned int)
{
  OutputDTITubeSpatialObjectPointer  outputDTITubeSpatialObject = OutputDTITubeSpatialObjectType::New();
  return dynamic_cast< DataObject *>( outputDTITubeSpatialObject.GetPointer() );
}




/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
void 
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::SetInput(unsigned int idx,const InputImageType *input)
{
  // process object is not const-correct, the const_cast
  // is required here.
  this->ProcessObject::SetNthInput(idx, 
                                   const_cast< InputImageType * >(input) );
}


  
/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
const typename ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>::InputImageType *
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::GetInput(unsigned int idx)  const
{
  return dynamic_cast<const InputImageType*>
    (this->ProcessObject::GetInput(idx));
}

/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
typename ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>::InputImageType *
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::GetInput(unsigned int idx) 
{
  return dynamic_cast<InputImageType*>
    (this->ProcessObject::GetInput(idx));
}

 
/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
const typename ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>::OutputDTITubeSpatialObjectType *
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::GetOutput(void)  const
{
  return dynamic_cast<const OutputDTITubeSpatialObjectType*>
    (this->ProcessObject::GetOutput(0));
}

/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
typename ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>::OutputDTITubeSpatialObjectType *
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::GetOutput(void) 
{
  return dynamic_cast<OutputDTITubeSpatialObjectType*>
    (this->ProcessObject::GetOutput(0));
}


/**
 *
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
void 
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}



/**
 * copy information from first input to all outputs
 * This is a void implementation to prevent the 
 * ProcessObject version to be called
 */
template <class TInputImage, class TOutputDTITubeSpatialObject>
void 
ImageToDTITubeSpatialObjectFilter<TInputImage,TOutputDTITubeSpatialObject>
::GenerateOutputInformation()
{
}


} // end namespace itk

#endif
