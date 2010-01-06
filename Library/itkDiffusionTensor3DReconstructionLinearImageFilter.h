/*=========================================================================

  Program:   NeuroLib (DTI command line tools)
  Language:  C++
  Date:      $Date: 2009/01/09 15:39:51 $
  Version:   $Revision: 1.3 $

  Author:    Casey Goodlett (gcasey@sci.utah.edu)

  Copyright (c)  Casey Goodlett. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDiffusionTensor3DReconstructionLinearImageFilter_h_
#define __itkDiffusionTensor3DReconstructionLinearImageFilter_h_

#include "itkDiffusionTensor3DReconstructionImageFilterBase.h"

namespace itk{
/** \class DiffusionTensor3DReconstructionLinearImageFilter
 * \brief This class derives from the
 * DiffusionTensor3DReconstructionImageFilterBase to implement a
 * least-squares tensor estimation.
 * 
 * \par References:
 * \li<a href="http://lmi.bwh.harvard.edu/papers/pdfs/2002/westinMEDIA02.pdf">[1]</a> 
 * <em>C.F.Westin, S.E.Maier, H.Mamata, A.Nabavi, F.A.Jolesz, R.Kikinis,
 * "Processing and visualization for Diffusion tensor MRI", Medical image
 * Analysis, 2002, pp 93-108.</em>
 * \li<a href="splweb.bwh.harvard.edu:8000/pages/papers/westin/ISMRM2002.pdf">[2]</a>
 * <em>A Dual Tensor Basis Solution to the Stejskal-Tanner Equations for DT-MRI</em>
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
 * DiffusionTensor3DReconstructionWeightedImageFilter
 * DiffusionTensor3DReconstructionRicianImageFilter
 * \ingroup Multithreaded  TensorObjects
 */
template< class TGradientImagePixelType,
          class TTensorPrecision=double >
class ITK_EXPORT DiffusionTensor3DReconstructionLinearImageFilter :
  public DiffusionTensor3DReconstructionImageFilterBase< TGradientImagePixelType,
                                                                TTensorPrecision>
{
public:
  typedef DiffusionTensor3DReconstructionLinearImageFilter       Self;
  typedef SmartPointer<Self>                                     Pointer;
  typedef SmartPointer<const Self>                               ConstPointer;
  typedef DiffusionTensor3DReconstructionImageFilterBase<TGradientImagePixelType,
    TTensorPrecision >                                           Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(Self, Superclass);

protected:
  DiffusionTensor3DReconstructionLinearImageFilter() {};
  virtual ~DiffusionTensor3DReconstructionLinearImageFilter() {};
  
  virtual vnl_vector< TTensorPrecision > 
    EstimateTensor(const vnl_vector<TTensorPrecision>& S) const;

};

}
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionLinearImageFilter.txx"
#endif

#endif
