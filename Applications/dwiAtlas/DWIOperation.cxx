//Some operations in this code are based on code from Casey Goodlett, Marc Niethammer, and Francois Budin

#include <iostream>
#include <fstream>
#include "DWIOperationCLP.h"
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

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

//custom filter
#include "SphericalHarmonicsInterpolationFilter.h"

//Legendre code from GNU Scientific Library
//#include <gsl/gsl_sf_legendre.h>

#define DIMENSION 3
#define PI 3.14159265358


//command-line parameters
struct params
{
    std::string dwiInputVolume;
    std::string dwiOutputVolume;
    std::string gradientVectorFile;
    int SHOrder;
    double lambda;
} ;


//figures out the type of each numerical value in the image
//This code from Francois Budin
void GetImageType (std::string fileName, itk::ImageIOBase::IOPixelType &pixelType, itk::ImageIOBase::IOComponentType &componentType)
{
    typedef itk::Image< unsigned char , 3 > ImageType ;
    itk::ImageFileReader< ImageType >::Pointer imageReader = itk::ImageFileReader< ImageType >::New();
    imageReader->SetFileName( fileName.c_str() ) ;
    imageReader->UpdateOutputInformation() ;
    pixelType = imageReader->GetImageIO()->GetPixelType() ;
    componentType = imageReader->GetImageIO()->GetComponentType() ;
}


//counts # of gradients in the provided text file.
int getNumOutputGradients(const std::string &gfname)
{
    int count = 0;

    std::ifstream inFile;
	inFile.open(gfname.c_str(), std::ios::in);

	if (!inFile)
	{
		std::cout << "Could not find gradient direction input file. Quitting.\n";
		exit(0);
	}
	else
	{
		std::string currentLine;
		
		while (!inFile.eof())
		{
			std::getline(inFile, currentLine);

			if ((currentLine.length() > 0) && !(currentLine.at(0) == '#')) //check if it's a comment line
			{
                count++;			
            }
        }
        inFile.close();
    }

    return count;
}


//fills up a vnl_matrix with the gradients in the provided text file
void parseGradientFile(const std::string &gfname, vnl_matrix<double> &newgradients)
{
    int count = 0;    
    int total = newgradients.rows();

    std::ifstream inFile;
	inFile.open(gfname.c_str(), std::ios::in);
	if (!inFile)
	{
		std::cout << "Could not find gradient direction input file. Quitting.\n";
		exit(0);
	}
	else
	{
		std::cout << "Reading new gradient directions...\n\n";
		
		std::string currentLine;
		
		while (!inFile.eof() && count < total)
		{
			std::getline(inFile, currentLine);

			if ((currentLine.length() > 0) && !(currentLine.at(0) == '#')) //check if it's a comment line
			{
                std::istringstream iss(currentLine);
                iss >> newgradients(count, 0) >> newgradients(count, 1) >> newgradients(count, 2);		
                count++;	
            }
        }
        inFile.close();
    }
}


//main processing method
template<class DWIPixelType> int process(params list)
{
    //we have a "vector image".  Can be thought of as a large 3-d volume, 
    //where every "pixel" is actually a _vector_ representing the n intensity values
    //for the n gradients, plus the baseline(s).    
    typedef itk::VectorImage<DWIPixelType,DIMENSION> DiffusionImageType;    

    //acquire information for input/output files
    const std::string infile = list.dwiInputVolume;
    const std::string outfile = list.dwiOutputVolume;

    //a reader for vectorimages ("DiffusionImages")
    typedef itk::ImageFileReader<DiffusionImageType> ImageFileReaderType;
    typename ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
    
    //make sure the infile works with the reader
    reader->SetFileName(infile);
    try
    {
        //update the reader to get output from the pipeline
        std::cout << "Reading input file...\n";
        reader->Update();
    }
    catch (itk::ExceptionObject & e)
    {
        std::cerr << "Something is wrong with your input file:\n";        
        std::cerr << e <<std::endl;
        return EXIT_FAILURE;
    }
    typename DiffusionImageType::Pointer dwimg = reader->GetOutput();

    //parse the new gradient directions in the input text file
    int numnewgrads = getNumOutputGradients(list.gradientVectorFile);
    vnl_matrix<double> newgrads(numnewgrads, 3, 0.0);
    parseGradientFile(list.gradientVectorFile, newgrads);

    //Filter the image
    typedef itk::SphericalHarmonicsInterpolationFilter<DiffusionImageType, DiffusionImageType> SphericalHarmonicsInterpolationFilterType; 
    typename SphericalHarmonicsInterpolationFilterType::Pointer shFilter = SphericalHarmonicsInterpolationFilterType::New();
    shFilter->SetInput(dwimg);
    shFilter->SetNewGradients(newgrads);
    shFilter->SetOrder(list.SHOrder);
    shFilter->SetLambda(list.lambda);
    shFilter->Update();

    int num_baselines = shFilter->GetNumBaselines();

    //********************    
    //METADATA DICTIONARY*
    //********************
    std::cout << "Constructing metadata for output file...\n";

    typedef itk::MetaDataDictionary DictionaryType;
    const DictionaryType & dictionary = reader->GetMetaDataDictionary();
  
    std::vector<std::string> imgMetaKeys = dictionary.GetKeys();
    std::vector<std::string>::iterator itKey = imgMetaKeys.begin();
    std::string metaString;

    // create a new metadata dictionary
    itk::MetaDataDictionary newDictionary;

    typedef itk::MetaDataObject< std::string > MetaDataStringType;
  
    DictionaryType::ConstIterator itr = dictionary.Begin();
    DictionaryType::ConstIterator end = dictionary.End();

    //copy over relevant info consistent between input/output images    
    while ( itr != end )
    {
        itk::MetaDataObjectBase::Pointer entry = itr->second;
        MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;

        if ( entryvalue )
        {
            std::string curKey = std::string(itr->first);
            std::string curVal = std::string(entryvalue->GetMetaDataObjectValue());

            if (!(curKey.substr(0, std::string("DWMRI_b-value").length()).compare(std::string("DWMRI_b-value")) == 0)
                && !(curKey.substr(0, std::string("DWMRI_gradient").length()).compare(std::string("DWMRI_gradient")) == 0))
            {
                itk::EncapsulateMetaData<std::string>(newDictionary, curKey, curVal);
            }
        }

        ++itr;
    }

    //write new gradient vectors (and b-value) into metadatadictionary
    std::string curKey = std::string("DWMRI_b-value");
    std::ostringstream bv;
    bv << shFilter->GetBValue();    
    std::string curVal = "" + bv.str();
    itk::EncapsulateMetaData<std::string>(newDictionary, curKey, curVal);

    char gnum[4]; //padded gradient direction id

    //take care of the baseline values
    for (int i = 0; i < num_baselines; i++)
    {
        sprintf(gnum, "%04d", i);
        curKey = "DWMRI_gradient_" + std::string(gnum);
        curVal = std::string("0 0 0");
        itk::EncapsulateMetaData<std::string>(newDictionary, curKey, curVal);
    }

    //take care of the rest of the gradient directions
    for (int i = 0; i < numnewgrads; i++)
    {      
        sprintf(gnum, "%04d", i+num_baselines);
        std::ostringstream g1;
        std::ostringstream g2;
        std::ostringstream g3;
        g1 << newgrads(i, 0);
        g2 << newgrads(i, 1);
        g3 << newgrads(i, 2);

        curKey = "DWMRI_gradient_" + std::string(gnum);
        curVal = " " + g1.str() + " " + g2.str() + " " + g3.str();

        itk::EncapsulateMetaData<std::string>(newDictionary, curKey, curVal);
    }
    //************************    
    //END METADATA DICTIONARY*
    //************************
    
    // write output
    itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
    io->SetFileTypeToBinary();    
    io->SetMetaDataDictionary(newDictionary);

    typedef itk::ImageFileWriter<DiffusionImageType> ImageFileWriterType;
    typename ImageFileWriterType::Pointer writer = ImageFileWriterType::New();

    writer->SetInput(shFilter->GetOutput());
    writer->UseInputMetaDataDictionaryOff();
    writer->SetImageIO(io);
  
    writer->SetFileName(outfile);
    writer->UseCompressionOn(); std::cout << "Output compression is turned ON.\n"; //comment/uncomment to enable compression
    try
    {
        std::cout << "Writing output to disk (this may take a few minutes)...\n";        
        writer->Update();
    }
    catch (itk::ExceptionObject & e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Finished.\n";

    return EXIT_SUCCESS;
}


int main(int argc, char* argv[])
{
    //takes care of generating usage message  
    PARSE_ARGS;

    //need to take infile, outfile at command-line
    //explain how to use program
    //Usage: argv[0] infile outfile inputgradients shorder
    //Description: Takes .nrrd image, reapproximates intensity values using
    //              spherical harmonics basis functions.    

    //take care of parameters    
    params list;
    list.dwiInputVolume = dwiInputVolume;
    list.dwiOutputVolume = dwiOutputVolume;
    list.gradientVectorFile = gradientVectorFile;
    list.SHOrder = SHOrder;
    list.lambda = lambda;

    //quick check for SHOrder
    if (list.SHOrder % 2 != 0)
    {
        std::cout << "The order you specified for the spherical harmonics (" << list.SHOrder
                  << ") needs to\n be an even number.\nExiting...\n";
        return EXIT_FAILURE;
    }

    //figure out what kind of pixels are used (i.e, RGB, RGBA, OFFSET, etc)
    //figure out whether representation is UCHAR, CHAR, ..., etc.)
    //see http://www.itk.org/Doxygen/html/classitk_1_1ImageIOBase.html
    //code from Francois Budin
    itk::ImageIOBase::IOPixelType pixelType ;
    itk::ImageIOBase::IOComponentType componentType ;
    GetImageType( list.dwiInputVolume, pixelType , componentType ) ;

    switch( componentType )
    {
        case itk::ImageIOBase::UCHAR:
            return process<unsigned char>(list);
            break;
        case itk::ImageIOBase::CHAR:
            return process<char>(list);
            break;
        case itk::ImageIOBase::USHORT:
            return process<unsigned short>(list);
            break;
        case itk::ImageIOBase::SHORT:
            return process<short>(list);
            break;
        case itk::ImageIOBase::UINT:
            return process<unsigned int>(list);
            break;
        case itk::ImageIOBase::INT:
            return process<int>(list);
            break;
        case itk::ImageIOBase::ULONG:
            return process<unsigned long>(list);
            break;
        case itk::ImageIOBase::LONG:
            return process<long>(list);
            break;
        case itk::ImageIOBase::FLOAT:
            return process<float>(list);
            break;
        case itk::ImageIOBase::DOUBLE:
            return process<double>(list);
            break;
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
            std::cout << "Unknown image component type.\nExiting...\n";
            break;
    }
    return 0;
}



