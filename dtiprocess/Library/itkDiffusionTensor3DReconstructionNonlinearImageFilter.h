/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.5 $

  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDiffusionTensor3DReconstructionNonlinearImageFilter_h_
#define __itkDiffusionTensor3DReconstructionNonlinearImageFilter_h_

#include "itkDiffusionTensor3DReconstructionImageFilterBase.h"
#include "vnl/vnl_least_squares_function.h"

namespace itk{

/** \class DiffusionTensor3DReconstructionNonlinearImageFilter
 * \brief This class derives from the
 * DiffusionTensor3DReconstructionImageFilterBase to implement a
 * non-linear least-squares tensor estimation.  This requires an
 * optimization and this implementation uses levenberg-marquardt.
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
 * DiffusionTensor3DReconstructionLinearImageFilter
 * DiffusionTensor3DReconstructionWeightedImageFilter
 * DiffusionTensor3DReconstructionRicianImageFilter
 * \ingroup Multithreaded  TensorObjects
 */
template< class TGradientImagePixelType,
          class TTensorPrecision=double >
class ITK_EXPORT DiffusionTensor3DReconstructionNonlinearImageFilter :
  public DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                                TTensorPrecision>
{
public:
  typedef DiffusionTensor3DReconstructionNonlinearImageFilter       Self;
  typedef SmartPointer<Self>                                     Pointer;
  typedef SmartPointer<const Self>                               ConstPointer;
  typedef DiffusionTensor3DReconstructionImageFilterBase<TGradientImagePixelType,
    TTensorPrecision >                                           Superclass;

  typedef typename Superclass::GradientPixelType                 GradientPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(Self, Superclass);

  /** Set step size for optimizer for non-linear fit.  A
  levenburg-Marquadt optimizer is used.  The default is 1.0e-10. */
  itkSetMacro( Step, double );
  /** Get step size for optimizer for non-linear fit.  A
  levenburg-Marquadt optimizer is used.  The default is 1.0e-10. */
  itkGetMacro( Step, double );

protected:
  DiffusionTensor3DReconstructionNonlinearImageFilter() : m_Step(1.0e-10) {};
  virtual ~DiffusionTensor3DReconstructionNonlinearImageFilter() {};
  
  virtual vnl_vector< TTensorPrecision > 
    EstimateTensor(const vnl_vector<TTensorPrecision>& S) const;
  
  double                                                 m_Step;
};

}
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionNonlinearImageFilter.txx"
#endif

#endif
