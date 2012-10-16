#include "sphericalHarmonicsFunctions.h"

#include <cmath>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <vector>

namespace sh
{

void initSH()
{
  createSHLookupMaps();
}

void createSHLookupMaps()
{
  // NEEDS TO BE CALLED BEFORE ANYTHING ELSE

  jtokmap.clear();
  jtommap.clear();

  int maxJ = (MAXK*MAXK + MAXK + 2)/2 + MAXK;

  jtokmap.resize( maxJ+1 ); // needs to be +1, because we index
			    // directly with the j-value
  jtommap.resize( maxJ+1 );

  for (int k = 0; k <= MAXK; k += 2)
    {
    for (int m = -k; m <= k; m++)
      {
      int jval = (k*k + k + 2)/2 + m;
      
      jtokmap[jval] = k;
      jtommap[jval] = m;

      }
    }
}

//returns the factorial of a number
// LEGACY: Should no longer be needed. Replaced by call to boost's
//factorial function
/*double
factorial(unsigned int n)
{
    if (n <= 1)
        return 1;
    else 
        return n * factorial(n-1);
}*/

//finds the maximum possible order for spherical harmonics that can be used
//note: order is even
int
findMaxOrder(unsigned int numgvectors)
{
    unsigned int numterms = 0;
    unsigned int order = 0;
    
    for (unsigned int ord = 2; ord <= 50; ord += 2)
    {
         numterms = (ord + 1)*(ord + 2) / 2;
        if (numterms > numgvectors)
        {
            break;
        }
        else
        {
            order = ord;
        }
    }
    
    return order;
}

//translate order of SH basis into number of terms
int
getNumTerms(unsigned int orderval)
{
    return ((orderval + 1)*(orderval + 2) / 2);
}

//convert cartesian coordinates to spherical coordinates (physics notation)

template <class RealType>
void
cart2sph(vnl_matrix<RealType> &gvectors, vnl_matrix<RealType> &sgvectors )
{
    RealType x, y, z, theta, phi, r, hypotxy;

    //based on matlab cart2sph m file
    for (unsigned int v = 0; v < gvectors.rows(); v++)
    {
        x = gvectors(v,0);
        y = gvectors(v,1);
        z = gvectors(v,2);

        //compute r
        hypotxy = sqrt(x*x + y*y);
        r = sqrt(x*x + y*y + z*z);

        //compute elev
        phi = atan2(z, hypotxy);

        //compute az 
        theta = atan2(y,x);

        //end up modifying and reversing names to fit our needs (physics notation)
        sgvectors(v,0) = (-1 * phi) + ((RealType)PI/2); //theta (physics)
        sgvectors(v,1) = theta; //phi (physics)
    }
}   

//return the appropriate spherical harmonic based on Legendre polynomials
//multiplication by sin or cos (for rpart or ipart) taken care of in calling function
template <class RealType>
RealType
getSphericalHarmonic(int l, int m, RealType theta, RealType phi)
{
  int m_abs = (int)fabs(m);    

  ////m_abs + 1 for matlab indexing (starting at 1)
  RealType Plmx = boost::math::legendre_p(l, m_abs, cos(theta)); //m_abs + 1 for matlab indexing (starting at 1)

  //double SR = sqrt(((2.0*l + 1.0)/4.0*PI)*((factorial(l -
  //m_abs))/(factorial(l + m_abs))));

  // TODO: make sure that this returns the correct values

  //RealType SR = sqrt(((2.0*l +
  //1.0)/4.0*(RealType)PI)*((boost::math::unchecked_factorial<int>(l-m_abs))/(boost::math::unchecked_factorial<int>(l+m_abs))));

  // TODO: Check if it is possible to use an int version of boost's
  // factorial function, seems like it crashes when I try to do this

  RealType SR = (RealType)(sqrt(((2.0*l + 1.0)/4.0*PI)*((boost::math::factorial<double>(l-m_abs))/(boost::math::factorial<double>(l+m_abs)))));
  
  RealType factor = SR * Plmx;
  
  return factor;     
}

//returns a value for a particular element of the SH basis matrix
template <class RealType>
RealType
getBasisMatrixValue(int curterm, RealType theta, RealType phi)
{
    //l is the 'order'
    //m is the 'phase factor'

    //let k = 0, 2, 4, ..., l (even number between zero and l inclusive)
    //then m = -k, ..., 0, ..., k

    //examples: (l = 4, in this case)
    //k = 0, m = 0; 
    //k = 2, m = {-2, -1, 0, 1, 2}; 
    //k = 4, m = {-4, -3, -2, -1, 0, 1, 2, 3, 4}

    //new indexing convention
    //j = (k^2 + k + 2) / 2

    //figure out the correct values of k and m for a specified j
    /*int kval = 0;
    int mval = 0;
    int jval = 0;

    for (int k = 0; k <= 20; k += 2)
    {
        for (int m = -k; m <= k; m++)
        {
            jval = (k*k + k + 2)/2 + m;

            if (jval == curterm)
            {
                kval = k;
                mval = m;

                break;
            }
        }
    }


    // TODO: replace this by a list and see if it speeds up things!!!

    kval = jtokmap[ curterm ];
    mval = jtommap[ curterm ];*/

    //set k and m to the appropriate values
  int k = jtokmap[ curterm ]; //kval;
  int m = jtommap[ curterm ]; //mval;

    //return the correct basis function value
    RealType Yj = 0;

    if ((-k <= m) && (m < 0))
    {
        //real part
        Yj = sqrt(2.0) * getSphericalHarmonic<RealType>(k, m, theta, phi) * cos(m*phi);
    }
    else if (m == 0)
    {
        Yj = getSphericalHarmonic<RealType>(k, 0, theta, phi);
    }
    else // (m > 0 && m <= k)
    {
        //imaginary part
        Yj = sqrt(2.0) * getSphericalHarmonic<RealType>(k, m, theta, phi) * sin(m*phi);
    }

    return Yj;
}

//generates the entire SH basis matrix
template <class RealType>
void
generateSHBasisMatrix(int numterms, int numgradients, vnl_matrix<RealType> &sgradients, vnl_matrix<RealType> &basismat)
{
    vnl_matrix<RealType> curgrad(1,2, 0.0);

    //basismat is the matrix we are modifying in-place.  It is numgradients x numterms in size.
    for (int i = 0; i < numgradients; i++) //each gradient vector
    {
        //extract current gradient
        curgrad(0, 0) = sgradients(i, 0); //theta
        curgrad(0, 1) = sgradients(i, 1); //phi
        
        for (int j = 1; j <= numterms; j++)
        {
            //fill in (i,j)th position in the basismat with the appropriate term
            RealType this_value = getBasisMatrixValue<RealType>(j, curgrad(0,0), curgrad(0,1));
            basismat(i,j-1) = this_value;
        }
    }
}

unsigned int getNumTermsAndCheckOrder( unsigned int order, unsigned int numoriggvectors, unsigned int & usedOrder )
{
  unsigned int max_order = findMaxOrder(numoriggvectors);
  
  //use max order if no order is provided at runtime (provide warning)
  //also use max order if provided order is too large.
  if (order <= 0)
    {    
    order = max_order - 2;
    std::cout << "WARNING: Using maxOrder-2 of Spherical harmonic basis\n functions, as no order was provided.\n";
    }
  else if (order > max_order)
    {
    order = max_order;
    std::cout << "WARNING: The provided order is too large; using maximum allowable order: " << order << ".\n";
    }
  else if (order == max_order)
    {
    std::cout << "Note: The order you provided (" << order << ") is the maximum allowable order.\n";
    }    

  usedOrder = order;
  
  //# of terms per row of the SH basis matrix    

  return getNumTerms(order);
}

template <class RealType>
void
computeSHRegularizationMatrix( vnl_diag_matrix<RealType> &regdiag, const unsigned int num_terms, const unsigned int order )
{
  // TODO: this could be optimized by just computing it once
  //std::cout << "Regularization term is being used. Lambda = " << lambda << std::endl;        

  //use the regularization term, if necessary
  vnl_vector<RealType> regvec(num_terms, 0.0);
    
  int currentOffset = 0;        
  for (unsigned int i = 0; i <= order; i+= 2)
    {
    int mult = 2*i + 1;
    int curVal = i;
    for (int j = currentOffset; j < (currentOffset+mult); j++)
      {
      regvec(j) = curVal;
      }
    currentOffset = currentOffset + mult;
    }
  
  vnl_vector<RealType> regvecsq(num_terms, 0.0);
  vnl_vector<RealType> regvecplus1sq(num_terms, 0.0);
  //square the vectors
  for (unsigned int i = 0; i < regvec.size(); i++)
    {
    regvecsq(i) = regvec(i) * regvec(i);
    regvecplus1sq = (regvec(i) + 1) * (regvec(i) + 1);
    }
  vnl_vector<RealType> regvecfinal(num_terms, 0.0);    
  //multiply
  for (unsigned int i = 0; i < regvec.size(); i++)
    {
    regvecfinal(i) = regvecsq(i) * regvecplus1sq(i);
    }
  
  //vnl_diag_matrix<RealType> regdiag(regvecfinal);
  regdiag.set(regvecfinal);
  
}

template <class RealType>
void 
generateSHBasisMatrixPseudoInverse( vnl_matrix<RealType>& cart_gradients, const RealType lambda, const unsigned int num_terms, const unsigned int order, vnl_matrix<RealType>& sh_basis_mat_pi )
{
  unsigned int numoriggvectors = cart_gradients.rows();

  //determine size of SH basis matrix
  vnl_matrix<RealType> sh_basis_mat(numoriggvectors, num_terms, 0.0);
  
  computeSHOrigBasisMat( sh_basis_mat, cart_gradients, num_terms );

  //determine size of SH basis mat pseudoinverse, and calculate it
  
  sh_basis_mat_pi.set_size(num_terms, num_terms);
  sh_basis_mat_pi = sh_basis_mat.transpose() * sh_basis_mat;
  
  if (lambda != 0)
    {

    vnl_diag_matrix<RealType> regdiag( num_terms );
    computeSHRegularizationMatrix( regdiag, num_terms, order );
    
    //add regularization term
    sh_basis_mat_pi = sh_basis_mat_pi + (regdiag.asMatrix() * lambda);
    }
  
  sh_basis_mat_pi = vnl_matrix_inverse<RealType>(sh_basis_mat_pi);
  sh_basis_mat_pi = sh_basis_mat_pi * sh_basis_mat.transpose(); 
  
}

template <class RealType>
RealType
huberWeightFcn( RealType scaledResidual, RealType C, bool &isOutlier )
{
  if ( fabs( scaledResidual ) <= C )
    {
    isOutlier = false;
    return 1;
    }
  else
    {
    isOutlier = true;
    return C/fabs(scaledResidual);
    }
}

template <class RealType>
void
computeSHHuberWeightMatrix( vnl_diag_matrix<RealType> &W, const unsigned int numoriggvectors, const vnl_vector<RealType>& shCoeffs, const vnl_matrix<RealType> &B, const vnl_vector<RealType> Y, const RealType C, const RealType sigma, unsigned int & nrOfOutliers )
{
  // IMPORTANT: Assumes that the signal lives in the log domain
  // estimation of the true signal is based on a first fit

  // compute the reconstructed values, given a set of spherical
  // harmonics coefficients

  vnl_vector<RealType> SE(numoriggvectors, 0.0);
  SE = B*shCoeffs;  // shCoeffs are the current spherical harmonics
		    // coefficients
  // now compute the residuals

  vnl_vector<RealType> residuals(numoriggvectors, 0.0 );
  residuals = SE-Y; // where Y are the measurements

  // now we can compute the weights based on the Huber function
  // (for very large values of C this becomes weighted least squares

  vnl_vector<RealType> wV(numoriggvectors, 0.0 );

  nrOfOutliers = 0;

  for ( unsigned int iI=0; iI<numoriggvectors; iI++ )
    {
    RealType SOrig = exp( SE(iI) );
    bool isOutlier;
    wV = SOrig*SOrig/(sigma*sigma)*huberWeightFcn(SOrig/sigma*residuals(iI), C, isOutlier );
    if ( isOutlier ) nrOfOutliers++;
    }
  
  W.set( wV );  // fill the entries of the weight matrix

}

template <class RealType>
void
computeSHOrigBasisMat( vnl_matrix<RealType>& sh_basis_mat, vnl_matrix<RealType>& cart_gradients, const unsigned int num_terms )
{
  unsigned int numoriggvectors = cart_gradients.rows();

  vnl_matrix<RealType> sph_gradients(numoriggvectors, 2);    
  //convert cartesian gradients to spherical coordinates    
  cart2sph<RealType>(cart_gradients, sph_gradients);
  
  //generate the SH basis matrix
  generateSHBasisMatrix<RealType>(num_terms, numoriggvectors, sph_gradients, sh_basis_mat);
}

template <class RealType>
void
generateWeightedSHBasisMatrixPseudoInversePrecomputed( const vnl_matrix<RealType> &B, const vnl_diag_matrix<RealType> &L, const vnl_diag_matrix<RealType>& W, const RealType lambda, vnl_matrix<RealType>& sh_basis_mat_pi )
{
  sh_basis_mat_pi = B.transpose() * W * B + L.asMatrix() * lambda;
  sh_basis_mat_pi = vnl_matrix_inverse<RealType>(sh_basis_mat_pi);
  sh_basis_mat_pi = sh_basis_mat_pi * B.transpose()*W; 
}

template <class RealType>
void
generateSHBasisMatrixPseudoInversePrecomputed( const vnl_matrix<RealType> &B, const vnl_diag_matrix<RealType> &L, const RealType lambda, vnl_matrix<RealType>& sh_basis_mat_pi )
{
  sh_basis_mat_pi = B.transpose() * B + L.asMatrix() * lambda;
  sh_basis_mat_pi = vnl_matrix_inverse<RealType>(sh_basis_mat_pi);
  sh_basis_mat_pi = sh_basis_mat_pi * B.transpose(); 
}

template <class RealType>
void 
generateWeightedSHBasisMatrixPseudoInverse( vnl_matrix<RealType>& cart_gradients, const vnl_diag_matrix<RealType>& W, const RealType lambda, const unsigned int num_terms, const unsigned int order, vnl_matrix<RealType>& sh_basis_mat_pi )
{
  unsigned int numoriggvectors = cart_gradients.rows();
  
  //determine size of SH basis matrix
  vnl_matrix<RealType> sh_basis_mat(numoriggvectors, num_terms, 0.0);
  
  computeSHOrigBasisMat( sh_basis_mat, cart_gradients, num_terms );
  
  //determine size of SH basis mat pseudoinverse, and calculate it
  
  sh_basis_mat_pi.set_size(num_terms, num_terms);
  sh_basis_mat_pi = sh_basis_mat.transpose() * W * sh_basis_mat;
  
  if (lambda != 0)
    {

    vnl_diag_matrix<RealType> regdiag( num_terms );
    computeSHRegularizationMatrix( regdiag, num_terms, order );
    
    //add regularization term
    sh_basis_mat_pi = sh_basis_mat_pi + (regdiag.asMatrix() * lambda);
    }
  
  sh_basis_mat_pi = vnl_matrix_inverse<RealType>(sh_basis_mat_pi);
  sh_basis_mat_pi = sh_basis_mat_pi * sh_basis_mat.transpose()*W; 
  
}

}
