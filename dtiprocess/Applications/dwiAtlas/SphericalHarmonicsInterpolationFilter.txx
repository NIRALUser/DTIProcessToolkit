/*=========================================================================



=========================================================================*/

#ifndef _SphericalHarmonicsInterpolationFilter_txx
#define _SphericalHarmonicsInterpolationFilter_txx
#include "SphericalHarmonicsInterpolationFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include <cmath>

//#include <tr1/cmath>

#include <boost/math/special_functions/legendre.hpp>

namespace itk
{

template <class TInputImage, class TOutputImage>
void 
SphericalHarmonicsInterpolationFilter<TInputImage, TOutputImage>
::GenerateOutputInformation() 
{
    Superclass::GenerateOutputInformation();

    this->GetOutput()->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
    this->GetOutput()->SetLargestPossibleRegion( this->GetInput()->GetLargestPossibleRegion() );
    //this->GetOutput()->SetNumberOfComponentsPerPixel( this->GetInput()->GetNumberOfComponentsPerPixel() );
    if (numPixelComponents == -1)
    {
        std::cout << "No set of new gradient vectors provided.\nExiting...\n";
        exit(0);
    }
    else
    {
        this->GetOutput()->SetNumberOfComponentsPerPixel(numPixelComponents);
    }
    this->GetOutput()->Allocate();

}


template <class TInputImage, class TOutputImage>
void 
SphericalHarmonicsInterpolationFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
    // call the superclass' implementation of this method
    Superclass::GenerateInputRequestedRegion();
  
    // get pointers to the input and output
    typename Superclass::InputImagePointer inputPtr = const_cast< TInputImage * >( this->GetInput() );
    typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
  
    if ( !inputPtr || !outputPtr )
    {
        return;
    }

    // get a copy of the input requested region (should equal the output
    // requested region)
    typename TInputImage::RegionType inputRequestedRegion;
    inputRequestedRegion = inputPtr->GetRequestedRegion();

    // crop the input requested region at the input's largest possible region
    if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
        inputPtr->SetRequestedRegion( inputRequestedRegion );
        return;
    }
    else
    {
        // Couldn't crop the region (requested region is outside the largest
        // possible region).  Throw an exception.

        // store what we tried to request (prior to trying to crop)
        inputPtr->SetRequestedRegion( inputRequestedRegion );
    
        // build an exception
        InvalidRequestedRegionError e(__FILE__, __LINE__);
        e.SetLocation(ITK_LOCATION);
        e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
        e.SetDataObject(inputPtr);
        throw e;
    }
}


//returns the factorial of a number
template< class TInputImage, class TOutputImage>
double
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::factorial(int n)
{
    if (n <= 1)
        return 1;
    else 
        return n * factorial(n-1);
}


//finds the maximum possible order for spherical harmonics that can be used
//note: order is even
template< class TInputImage, class TOutputImage>
int
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::findMaxOrder(int numgvectors)
{
    int numterms = 0;
    int order = 0;
    
    for (int ord = 2; ord <= 50; ord += 2)
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
template< class TInputImage, class TOutputImage>
int
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::getNumTerms(int orderval)
{
    return ((orderval + 1)*(orderval + 2) / 2);
}


//convert cartesian coordinates to spherical coordinates (physics notation)
template< class TInputImage, class TOutputImage>
void
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::cart2sph(vnl_matrix<double> &gvectors, vnl_matrix<double> &sgvectors )
{
    double x, y, z, theta, phi, r, hypotxy;

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
        sgvectors(v,0) = (-1 * phi) + (PI/2); //theta (physics)
        sgvectors(v,1) = theta; //phi (physics)
    }
}   


//return the appropriate spherical harmonic based on Legendre polynomials
//multiplication by sin or cos (for rpart or ipart) taken care of in calling function
template< class TInputImage, class TOutputImage>
double
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::getSphericalHarmonic(int l, int m, double theta, double phi)
{
  int m_abs = (int)fabs(m);    

  //double Plmx = gsl_sf_legendre_Plm(l, m_abs, cos(theta)); //m_abs + 1 for matlab indexing (starting at 1)
  
  //double Plmx = std::tr1::assoc_legendre(l, m_abs, cos(theta));
  ////m_abs + 1 for matlab indexing (starting at 1)
  double Plmx = boost::math::legendre_p(l, m_abs, cos(theta)); //m_abs + 1 for matlab indexing (starting at 1)

  double SR = sqrt(((2.0*l + 1.0)/4.0*PI)*((factorial(l - m_abs))/(factorial(l + m_abs))));
  
  double factor = SR * Plmx;
  
  return factor;     
}


//returns a value for a particular element of the SH basis matrix
template< class TInputImage, class TOutputImage>
double
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::getBasisMatrixValue(int curterm, double theta, double phi)
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
    int kval = 0;
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

    //set k and m to the appropriate values
    int k = kval;
    int m = mval;

    //return the correct basis function value
    double Yj = 0;

    if ((-k <= m) && (m < 0))
    {
        //real part
        Yj = sqrt(2.0) * getSphericalHarmonic(k, m, theta, phi) * cos(m*phi);
    }
    else if (m == 0)
    {
        Yj = getSphericalHarmonic(k, 0, theta, phi);
    }
    else // (m > 0 && m <= k)
    {
        //imaginary part
        Yj = sqrt(2.0) * getSphericalHarmonic(k, m, theta, phi) * sin(m*phi);
    }

    return Yj;
}


//generates the entire SH basis matrix
template< class TInputImage, class TOutputImage>
void
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::generateSHBasisMatrix(int numterms, int numgradients, vnl_matrix<double> &sgradients, vnl_matrix<double> &basismat)
{
    vnl_matrix<double> curgrad(1,2, 0.0);

    //basismat is the matrix we are modifying in-place.  It is numgradients x numterms in size.
    for (int i = 0; i < numgradients; i++) //each gradient vector
    {
        //extract current gradient
        curgrad(0, 0) = sgradients(i, 0); //theta
        curgrad(0, 1) = sgradients(i, 1); //phi
        
        for (int j = 1; j <= numterms; j++)
        {
            //fill in (i,j)th position in the basismat with the appropriate term
            double this_value = getBasisMatrixValue(j, curgrad(0,0), curgrad(0,1));
            basismat(i,j-1) = this_value;
        }
    }
}

template< class TInputImage, class TOutputImage>
int
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage>
::countBaselines()
{
    int numbaselines = 0;    

    typename InputImageType::ConstPointer inputimg  = this->GetInput();

    //get access to the metadata dictionary (reference to the metadatadictionary)
    itk::MetaDataDictionary inputdict = inputimg->GetMetaDataDictionary();

    //get all the keys from the dictionary (ids for the attributes)
    std::vector<std::string> keys = inputdict.GetKeys();

    //get relevant information out of the metadata dictionary    
    for(std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it)
    {
        //count baselines
        if( it->find("DWMRI_gradient") != std::string::npos) //we found a gradient vector
        {
            std::string value;
            itk::ExposeMetaData<std::string>(inputdict, *it, value);
            std::istringstream iss(value);
            vnl_vector_fixed<double, 3> g;
            iss >> g[0] >> g[1] >> g[2];
            if (g[0] == 0 && g[1] == 0 && g[2] == 0)
                numbaselines++;
        }
    }
    std::cout << "Method determined " << numbaselines << " baselines.\n";
    return numbaselines;
}


template< class TInputImage, class TOutputImage>
void
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::GenerateData()
{
    // Get the input and output images
    typename OutputImageType::Pointer outputimg = this->GetOutput();
    typename InputImageType::ConstPointer inputimg  = this->GetInput();

    //**************************************************************
    // Get original gradient vectors out of the input image metadata
    //
    //**************************************************************    
    // Parse gradient directions from image header as specified by the
    // namic conventions defined at http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
    typedef vnl_vector_fixed<double, 3> GradientType;  
    typedef itk::VectorContainer<unsigned int, GradientType> GradientListType;
    typename GradientListType::Pointer gradientContainer = GradientListType::New();
    
    //get access to the metadata dictionary
    itk::MetaDataDictionary inputdict = inputimg->GetMetaDataDictionary();
    std::vector<std::string> keys = inputdict.GetKeys();   
    for(std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); ++it)
    {
        //obtain the b-value        
        if (it->find("DWMRI_b-value") != std::string::npos) //we found the b value
        {
            std::string bstr; 
            itk::ExposeMetaData<std::string>(inputdict, *it, bstr);
            std::istringstream iss(bstr); 
            iss >> bvalue; 
        }
        
        //obtain all the gradient vectors
        if( it->find("DWMRI_gradient") != std::string::npos) //we found a gradient vector
        {
            std::string value;
            itk::ExposeMetaData<std::string>(inputdict, *it, value);
            std::istringstream iss(value);
            GradientType g;
            iss >> g[0] >> g[1] >> g[2];
            unsigned int ind;
            std::string temp = it->substr(it->find_last_of('_')+1);
            ind = atoi(temp.c_str());
            //ind is the current gradient vector index
            gradientContainer->InsertElement(ind,g); 
        }
    }
    //*****************************************************************
    // End acquisition of original gradient vectors
    //
    //*****************************************************************

    int num_detected_vectors = inputimg->GetVectorLength();
    int baselinecounter = 0;
    for (int v = 0; v < num_detected_vectors; v++)
    {
        GradientType curgrad = gradientContainer->GetElement(v);
        if (curgrad[0] == 0 && curgrad[1] == 0 && curgrad[2] == 0)
        {
            baseline_indices.insert(v);
            baselinecounter++;
        }
        
    }

    int numoriggvectors = num_detected_vectors - baselinecounter;
    num_baselines = baselinecounter;

    //Set number of components per pixel in output
    //numPixelComponents = newgvectors.rows() + num_baselines;

    //create matrices for gradients
    vnl_matrix<double> cart_gradients(numoriggvectors, 3);
    vnl_matrix<double> sph_gradients(numoriggvectors, 2);    

    //fill up the cart_gradients with real values
    int insertioncounter = 0;
    for (int v = 0; v < num_detected_vectors; v++)
    {
        if (baseline_indices.count(v) == 0) //actual gradient, not baseline
        {        
            GradientType curgrad = gradientContainer->GetElement(v);
            cart_gradients(insertioncounter,0) = curgrad[0];
            cart_gradients(insertioncounter,1) = curgrad[1];
            cart_gradients(insertioncounter,2) = curgrad[2];
            insertioncounter++;
        }        
    } 
    
    //verify minimum number of gradient directions is present
    if (baselinecounter < 1 || numoriggvectors < 6)
    {
        std::cerr << "The input volume does not have enough gradient directions.\n";
        std::cerr << "Apart from at least 1 baseline, 6 additional directions are needed\n";
        std::cerr << "  for a total of 7 (minimum).\n";
        std::cerr << "Your input volume has only " << numoriggvectors << " non-baseline directions, and " << baselinecounter << "baselines.\n";
        exit(0); //TODO
    }
    
    //*********************************************************
    // Use original gradients to calculate pseudoinverse matrix
    //
    //*********************************************************    
    
    int max_order = findMaxOrder(numoriggvectors);

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

    //# of terms per row of the SH basis matrix    
    int num_terms = getNumTerms(order);

    std::cout << numoriggvectors << " gradient directions provided w/input image (not including baselines):\n";
    std::cout << cart_gradients << std::endl;
    std::cout << newgvectors.rows() << " new gradient directions provided for output (not including baselines):\n";
    std::cout << newgvectors << std::endl;
    std::cout << "b-value: " << bvalue << std::endl;
    std::cout << "Order for spherical harmonics interpolation: " << order << std::endl << std::endl;

    //convert cartesian gradients to spherical coordinates    
    cart2sph(cart_gradients, sph_gradients);

    //determine size of SH basis matrix
    vnl_matrix<double> sh_basis_mat(numoriggvectors, num_terms, 0.0);

    //generate the SH basis matrix
    generateSHBasisMatrix(num_terms, numoriggvectors, sph_gradients, sh_basis_mat);

    //determine size of SH basis mat pseudoinverse, and calculate it
    
    //use the regularization term, if necessary
    vnl_vector<double> regvec(num_terms, 0.0);

    vnl_matrix<double> sh_basis_mat_pi(num_terms, num_terms, 0.0);
    sh_basis_mat_pi = sh_basis_mat.transpose() * sh_basis_mat;

    if (lambda != 0)
    {
        std::cout << "Regularization term is being used. Lambda = " << lambda << std::endl;        
        
        int currentOffset = 0;        
        for (int i = 0; i <= order; i+= 2)
        {
            int mult = 2*i + 1;
            int curVal = i;
            for (int j = currentOffset; j < (currentOffset+mult); j++)
            {
                regvec(j) = curVal;
            }
            currentOffset = currentOffset + mult;
        }
        
        vnl_vector<double> regvecsq(num_terms, 0.0);
        vnl_vector<double> regvecplus1sq(num_terms, 0.0);
        //square the vectors
        for (unsigned int i = 0; i < regvec.size(); i++)
        {
            regvecsq(i) = regvec(i) * regvec(i);
            regvecplus1sq = (regvec(i) + 1) * (regvec(i) + 1);
        }
        vnl_vector<double> regvecfinal(num_terms, 0.0);    
        //multiply
        for (unsigned int i = 0; i < regvec.size(); i++)
        {
            regvecfinal(i) = regvecsq(i) * regvecplus1sq(i);
        }

        vnl_diag_matrix<double> regdiag(regvecfinal);
        
        //add regularization term
        sh_basis_mat_pi = sh_basis_mat_pi + (regdiag.asMatrix() * lambda);
    }
        
    sh_basis_mat_pi = vnl_matrix_inverse<double>(sh_basis_mat_pi);
    sh_basis_mat_pi = sh_basis_mat_pi * sh_basis_mat.transpose();   
    //************************************************************
    // End pseudoinverse matrix calculation
    //
    //************************************************************ 

    //create iterators to iterate over the vector elements of the 3-d images
    typedef itk::ImageRegionConstIterator< TInputImage > ConstIteratorType;
    ConstIteratorType iit(inputimg, inputimg->GetLargestPossibleRegion());
    typedef itk::ImageRegionIterator<TOutputImage> IteratorType;
    IteratorType oit(outputimg, outputimg->GetLargestPossibleRegion());

    //convert newgrads to spherical coordinates
    int numnewgvectors = newgvectors.rows();
    vnl_matrix<double> newgvectorssph(numnewgvectors, 2);
    cart2sph(newgvectors, newgvectorssph);

    //construct SHBnew
    vnl_matrix<double> sh_basis_mat_new(numnewgvectors, num_terms, 0.0);
    generateSHBasisMatrix(num_terms, numnewgvectors, newgvectorssph, sh_basis_mat_new);

    //create a matrix C for the coefficients
    vnl_matrix<double> C(num_terms, 1, 0.0);
    
    //create matrix to hold s_values (input image)    
    vnl_matrix<double> s_values(numoriggvectors, 1, 0.0);
    
    //create matrix to hold s_values (output image)
    vnl_matrix<double> s_values_new(numnewgvectors, 1, 0.0);
    
    //initialize a place to keep the baseline value at each element
    std::vector<ImageValueType> baselines;
    for (int i = 0; i < num_baselines; i++)
        baselines.push_back(0);
    
    //length of a vector element in the output image   
    int outputimg_vectorlength = numnewgvectors + baselinecounter;
    outputimg->SetVectorLength(outputimg_vectorlength); //TODO: does this actually do anything useful? 

    int blinsertioncounter = 0;

    std::cout << "Constructing output image...\n";  

    //input and output images are the same size, vector elements might be different sizes    
    for (iit.GoToBegin(), oit.GoToBegin(); !iit.IsAtEnd() && !oit.IsAtEnd(); ++iit, ++oit)    
    {
        //get current vector element from input image        
        InputImagePixelType this_in_element = iit.Get();     
        
        //create new vector element for output image (size it correctly)
        OutputImagePixelType this_out_element(outputimg_vectorlength);        

        //store baselines and add input s values into s values vector
        insertioncounter = 0;
        blinsertioncounter = 0;
        for (int v = 0; v < num_detected_vectors; v++)
        {
            if (baseline_indices.count(v) > 0) //baseline value, not gradient     
            {
                baselines[blinsertioncounter] = this_in_element[v];
                blinsertioncounter++;
            }
            else
            {
                s_values(insertioncounter, 0) = this_in_element[v];
                insertioncounter++;
            }     
        }

        //solve for C
        C = sh_basis_mat_pi * s_values;

        //generate intensity values
        s_values_new = sh_basis_mat_new * C;

        //put the values in the output image (baselines go in first positions)
        for (unsigned int v = 0; v < baselines.size(); v++)
        {
            this_out_element[v] = baselines[v];
        }
        insertioncounter = 0;        
        for (int v = baselines.size(); v < outputimg_vectorlength; v++)
        {
            double current_val = s_values_new(insertioncounter,0);           
            this_out_element[v] = (ImageValueType) current_val;
            insertioncounter++;
        }
        
        //place output vector element into output image
        oit.Set(this_out_element);
    }

    std::cout << "Output image construction complete.\n";
}


/**
 * Standard "PrintSelf" method
 */
template< class TInputImage, class TOutputImage>
void
SphericalHarmonicsInterpolationFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk


#endif

