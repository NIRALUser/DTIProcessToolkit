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
#include <vnl/algo/vnl_levenberg_marquardt.h>

namespace itk
{

namespace Functor
{
template<class TTensorPrecision>
class ExponentialLeastSquaresFunction : public vnl_least_squares_function
{
public:
  ExponentialLeastSquaresFunction(unsigned int                   number_of_residuals,
                                  const vnl_matrix<TTensorPrecision>* design_matrix,
                                  const vnl_vector<TTensorPrecision>* signal) :
    vnl_least_squares_function(7,number_of_residuals,no_gradient), 
    m_Design_Matrix(design_matrix),
    m_Signal(signal)
  {
  }
  
  void f(const vnl_vector<double> &x,
         vnl_vector<double> &fx) 
  {
    vnl_vector<double> predictedsignal(n_);
    predictedsignal = (*m_Design_Matrix) * x;
    for(unsigned int i = 0; i < n_; ++i)
    {    
      predictedsignal[i] = exp(predictedsignal[i]);
      fx[i] = (*m_Signal)[i] - predictedsignal[i];
    }
  }

private:
  const vnl_matrix<TTensorPrecision>* m_Design_Matrix;
  const vnl_vector<TTensorPrecision>* m_Signal;
};
}

template< class TGradientImagePixelType, class TTensorPrecision >
vnl_vector<TTensorPrecision>
DiffusionTensor3DReconstructionNonlinearImageFilter< TGradientImagePixelType, 
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

  vnl_vector< TTensorPrecision > estimate(this->m_TensorBasis * B);

  typedef Functor::ExponentialLeastSquaresFunction<TTensorPrecision> FittingFunctionType;
  
  unsigned int ng = this->m_NumberOfGradientDirections;
  FittingFunctionType lsqf(ng, 
                           &this->m_BMatrix, 
                           &S);
  
  // Levenberg-Marquardt is optimal for functions
  // which are sum of squared residuals
  vnl_levenberg_marquardt optimizer(lsqf);
  optimizer.set_x_tolerance(m_Step);
  optimizer.set_f_tolerance(1.0e-6);

  std::cout << this->m_BMatrix.size() << std::endl;
  if(!optimizer.minimize(estimate))
  {
    throw itk::ExceptionObject("Error in non-linear tensor estimation");
  }

  return estimate;
}

}
