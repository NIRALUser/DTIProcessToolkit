/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionWeightedImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionWeightedImageFilter_h_
#define __itkDiffusionTensor3DReconstructionWeightedImageFilter_h_

#include "itkDiffusionTensor3DReconstructionImageFilterBase.h"

namespace itk{

/** \class DiffusionTensor3DReconstructionLinearImageFilter      \
 *
 * \brief
 * This class derives from the
 * DiffusionTensor3DReconstructionImageFilterBase to implement a
 * weighted least-squares tensor estimation.
 * 
 * \par References:
 *
 * <em> Raymond Salvador and Alonso Pena and David K.  Menon and
                   T. Adrian Carpenter and John D. Pickard and Ed
                   T. Bullmore, "Formal characterization and extension of the
                   linearized diffusion tensor model", Human Brain
                   Mapping 24(2), 2005, pp 144-155 </em>
 * 
 * \note
 * This work is part of the National Alliance for Medical image Computing 
 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
 * for Medical Research, Grant U54 EB005149.
 *
 * \author Casey Goodlett.  Thanks to Xiaodong Tao, GE, for contributing the original
 * version of this class.
 * 
 * \sa DiffusionTensor3D SymmetricSecondRankTensor
 * DiffusionTensor3DReconstructionImageFilterBase
 * DiffusionTensor3DReconstructionNonlinearImageFilter
 * DiffusionTensor3DReconstructionLinearImageFilter
 * DiffusionTensor3DReconstructionRicianImageFilter
 * \ingroup Multithreaded  TensorObjects
 */
template< class TGradientImagePixelType, 
          class TTensorPrecision=double >
class ITK_EXPORT DiffusionTensor3DReconstructionWeightedImageFilter :
  public DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                         TTensorPrecision>
{

public:

  typedef DiffusionTensor3DReconstructionWeightedImageFilter Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef DiffusionTensor3DReconstructionImageFilterBase<TGradientImagePixelType,
    TTensorPrecision >                                           Superclass;
  
   /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(DiffusionTensor3DReconstructionWeightedImageFilter, 
                                                      ImageToImageFilter);
 
  /** Number of reweighting iterations to use.  The default and
   * recommmended number is 1.  */
  itkSetMacro( NumberOfIterations, unsigned int );
  itkGetMacro( NumberOfIterations, unsigned int );
  
protected:
  DiffusionTensor3DReconstructionWeightedImageFilter();
  virtual ~DiffusionTensor3DReconstructionWeightedImageFilter() {};

  virtual vnl_vector< TTensorPrecision > 
    EstimateTensor(const vnl_vector<TTensorPrecision>& S) const;

private:
  /** Number of reweighting iterations to use.  The default and
   * recommmended number is 1.  */
  unsigned int                                 m_NumberOfIterations;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionWeightedImageFilter.txx"
#endif

#endif

