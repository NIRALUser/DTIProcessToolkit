/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorBSplineInterpolateImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008/07/02 15:54:54 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorBSplineInterpolateImageFunction_txx
#define __itkVectorBSplineInterpolateImageFunction_txx

#include "itkVectorBSplineInterpolateImageFunction.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep, class TCoefficientType>
const unsigned long
VectorBSplineInterpolateImageFunction< TInputImage, TCoordRep , TCoefficientType >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep, class TCoefficientType>
VectorBSplineInterpolateImageFunction< TInputImage, TCoordRep , TCoefficientType >
::VectorBSplineInterpolateImageFunction()
{
  m_ComponentInterpolators.resize(Dimension);
  m_ComponentAdaptors.resize(Dimension);

  for ( unsigned int i = 0; i < Dimension; i++)
    {
    m_ComponentInterpolators[i] = ComponentInterpolateFunctionType::New();
    m_ComponentAdaptors[i] = ComponentAdaptorType::New();
    }
}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep, class TCoefficientType>
void
VectorBSplineInterpolateImageFunction< TInputImage, TCoordRep , TCoefficientType >
::PrintSelf(std::ostream& os, Indent indent) const
{
  os << "Vector BSpline" << std::endl;
  this->Superclass::PrintSelf(os,indent);
}

template <class TImageType, class TCoordRep, class TCoefficientType>
void 
VectorBSplineInterpolateImageFunction<TImageType,TCoordRep,TCoefficientType>
::SetInputImage(const TImageType * inputData)
{
  // Call super class input set
  this->VectorInterpolateImageFunction<TImageType,TCoordRep>::SetInputImage(inputData);

  if(inputData)
    {
    // data is now stored in m_Image
    for( unsigned int i = 0; i < Dimension; i++)
      {
      m_ComponentAdaptors[i]->SetInput(inputData);
      m_ComponentAdaptors[i]->SetIndex(i);
      m_ComponentAdaptors[i]->Update();
      
      m_ComponentInterpolators[i]->SetInputImage(m_ComponentAdaptors[i]->GetOutput());
      }
    }
}

/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep, class TCoefficientType>
typename VectorBSplineInterpolateImageFunction< TInputImage, TCoordRep , TCoefficientType >
::OutputType
VectorBSplineInterpolateImageFunction< TInputImage, TCoordRep , TCoefficientType >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType& index) const
{
  OutputType output;
  output.Fill(0.0);

  for( unsigned int i = 0; i < Dimension; i++ )
    {
    output[i] = m_ComponentInterpolators[i]->EvaluateAtContinuousIndex(index);
    }
  return output;
}

} // end namespace itk

#endif
