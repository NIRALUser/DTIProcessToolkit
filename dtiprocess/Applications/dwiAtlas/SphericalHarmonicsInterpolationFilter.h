#ifndef __SphericalHarmonicsInterpolationFilter_h
#define __SphericalHarmonicsInterpolationFilter_h

//this code builds on the ImageToImageFilter.
#include "itkImageToImageFilter.h"
#include "itkImage.h"

//Legendre code from GNU Scientific Library
//#include <gsl/gsl_sf_legendre.h>

#include <string.h>
#include <math.h>
#include <vector>
#include <set>

//itk includes go here
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include <itkNthElementImageAdaptor.h>
#include <itkImageToVectorImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVersion.h>

#define DIMENSION 3
#define PI 3.14159265358

namespace itk
{
/** \class SphericalHarmonicsInterpolationFilter
 *
 * \sa Image
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT SphericalHarmonicsInterpolationFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Extract dimension from input image. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage InputImageType;
    typedef typename    InputImageType::Pointer    InputImagePointer;
    typedef typename    InputImageType::RegionType InputImageRegionType; 
    typedef typename    InputImageType::PixelType  InputImagePixelType; //variable sized vector of valuetype
    typedef typename    InputImageType::PixelType::ValueType ImageValueType; //USE THIS TO DETERMINE INDIVIDUAL VECTOR VALUES
    typedef TOutputImage OutputImageType;
    typedef typename     OutputImageType::Pointer    OutputImagePointer;
    typedef typename     OutputImageType::RegionType OutputImageRegionType;
    typedef typename     OutputImageType::PixelType  OutputImagePixelType;

    /** Standard class typedefs. */
    typedef SphericalHarmonicsInterpolationFilter Self;    
    typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self>  ConstPointer;

    /** Typedef to describe the output and input image index and size types. */
    typedef typename TOutputImage::IndexType OutputImageIndexType;
    typedef typename TInputImage::IndexType  InputImageIndexType;
    typedef typename TOutputImage::SizeType  OutputImageSizeType;
    typedef typename TInputImage::SizeType   InputImageSizeType;
    typedef InputImageSizeType SizeType;
  
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SphericalHarmonicsInterpolationFilter, ImageToImageFilter);

    //TODO:    
    /** Extra public functions go here */
    void SetOrder(int neworder) { order = neworder; };
    void SetLambda(double lmbda) { lambda = lmbda; };
    void SetNewGradients(vnl_matrix<double> &newgrads) { newgvectors = newgrads; numPixelComponents = (newgrads.rows() + countBaselines()); };
    //void SetNewGradients(vnl_matrix<double> &newgrads) { newgvectors = newgrads; };
    int GetBValue() { return bvalue; };
    int GetNumBaselines() { return num_baselines; };


#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro(SameDimensionCheck,
        (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
    itkConceptMacro(InputConvertibleToOutputCheck,
        (Concept::Convertible<InputImagePixelType, OutputImagePixelType>));
  /** End concept checking */
#endif


protected:
    SphericalHarmonicsInterpolationFilter()
    {
      //newgvectors = NULL;
      newgvectors.clear();
        order = -1;
        bvalue = -1;
        numPixelComponents = -1;
        lambda = 0;
    }
    virtual ~SphericalHarmonicsInterpolationFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void GenerateOutputInformation();
    void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);
    void GenerateData();

private:
    SphericalHarmonicsInterpolationFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    //extra private functions go here
    double factorial(int n);
    int findMaxOrder(int numgvectors);
    int getNumTerms(int orderval);
    void cart2sph(vnl_matrix<double> &gvectors, vnl_matrix<double> &sgvectors );
    double getSphericalHarmonic(int l, int m, double theta, double phi);
    double getBasisMatrixValue(int curterm, double theta, double phi);
    void generateSHBasisMatrix(int numterms, int numgradients, vnl_matrix<double> &sgradients, vnl_matrix<double> &basismat);
    int countBaselines();

    //extra private variables go here
    vnl_matrix<double> newgvectors;
    int order;
    double lambda;
    int bvalue;
    int numPixelComponents;
    int num_baselines;
    std::set<int> baseline_indices;

};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SphericalHarmonicsInterpolationFilter.txx"
#endif

#endif

