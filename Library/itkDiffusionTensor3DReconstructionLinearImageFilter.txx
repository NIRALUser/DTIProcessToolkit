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
namespace itk
{

template< class TGradientImagePixelType, class TTensorPrecision >
vnl_vector<TTensorPrecision>
DiffusionTensor3DReconstructionLinearImageFilter< TGradientImagePixelType, 
                                                  TTensorPrecision >
::EstimateTensor(const vnl_vector<TTensorPrecision>& S) const
{
  vnl_vector< TTensorPrecision > B(this->m_NumberOfGradientDirections);
  for(unsigned int i = 0; i < S.size(); ++i)
  {
    if(S[i] == 0)
      B[i] = 0;
    else
      B[i] = log(S[i]);
  }

  return this->m_TensorBasis * B;
}

}
