#ifndef _SPHERICAL_HARMONICS_FUNCTIONS_
#define _SPHERICAL_HARMONICS_FUNCTIONS_

#include <itkArray2D.h>
#include <itkArray.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_diag_matrix.h>
#include <vector>

namespace sh
{

const double PI = 3.14159265358979323846;

std::vector<int> jtokmap;
std::vector<int> jtommap;

const int MAXK = 50;

//double factorial( unsigned int n );
int findMaxOrder( unsigned int numgvectors );
int getNumTerms(unsigned int orderval);

template <class RealType>
void cart2sph(vnl_matrix<RealType> &gvectors, vnl_matrix<RealType> &sgvectors );

template <class RealType>
RealType getSphericalHarmonic(int l, int m, RealType theta, RealType phi);

template <class RealType>
RealType getBasisMatrixValue(int curterm, RealType theta, RealType phi);

template <class RealType>
void generateSHBasisMatrix(int numterms, int numgradients, vnl_matrix<RealType> &sgradients, vnl_matrix<RealType> &basismat);

template <class RealType>
void generateSHBasisMatrixPseudoInverse( vnl_matrix<RealType>& cart_gradients, const RealType lambda, const unsigned int num_terms, const unsigned int order, vnl_matrix<RealType>& sh_basis_mat_pi );

template <class RealType>
void 
generateWeightedSHBasisMatrixPseudoInverse( vnl_matrix<RealType>& cart_gradients, const vnl_diag_matrix<RealType>& W, const RealType lambda, const unsigned int num_terms, const unsigned int order, vnl_matrix<RealType>& sh_basis_mat_pi );

template <class RealType>
void
generateWeightedSHBasisMatrixPseudoInversePrecomputed( const vnl_matrix<RealType> &B, const vnl_diag_matrix<RealType> &L, const vnl_diag_matrix<RealType>& W, const RealType lambda, vnl_matrix<RealType>& sh_basis_mat_pi );

template <class RealType>
void
generateSHBasisMatrixPseudoInversePrecomputed( const vnl_matrix<RealType> &B, const vnl_diag_matrix<RealType> &L, const RealType lambda, vnl_matrix<RealType>& sh_basis_mat_pi );

template <class RealType>
void
computeSHOrigBasisMat( vnl_matrix<RealType>& sh_basis_mat, vnl_matrix<RealType>& cart_gradients, const unsigned int num_terms );

template <class RealType>
void
computeSHRegularizationMatrix( vnl_diag_matrix<RealType> &regdiag, const unsigned int num_terms, const unsigned int order );

template <class RealType>
void
computeSHHuberWeightMatrix( vnl_diag_matrix<RealType> &W, const unsigned int numoriggvectors, const vnl_vector<RealType>& shCoeffs, const vnl_matrix<RealType> &B, const vnl_vector<RealType> Y, const RealType C, const RealType sigma, unsigned int &nrOfOutliers );

template <class RealType>
RealType
huberWeightFcn( RealType scaledResidual, RealType C, bool &isOutlier );


unsigned int getNumTermsAndCheckOrder( unsigned int order, unsigned int numoriggvectors, unsigned int & usedOrder );

void initSH();
void createSHLookupMaps();

}

#include "sphericalHarmonicsFunctions.txx"

#endif
