#include "auxFunctions.h"

#include <iostream>
#include <fstream>
#include <sstream>

// contains a set of auxiliary functions

namespace aux
{

template <class RealType>
void ConstructOutputMetaDataDictionary( itk::MetaDataDictionary & newDictionary, const itk::MetaDataDictionary & dictionary, unsigned int num_baselines, const vnl_matrix<RealType> &newgrads  )
{

  unsigned int num_newgrads = newgrads.rows();

  typedef itk::MetaDataDictionary DictionaryType;
  
  std::vector<std::string> imgMetaKeys = dictionary.GetKeys();
  std::vector<std::string>::iterator itKey = imgMetaKeys.begin();
  std::string metaString;
  
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
      
      if ( !(curKey.substr(0, std::string("DWMRI_gradient").length()).compare(std::string("DWMRI_gradient")) == 0))
	{
	itk::EncapsulateMetaData<std::string>(newDictionary, curKey, curVal);
	}
      }
    
    ++itr;
    }
  
  
  //write new gradient vectors (and b-value) into metadatadictionary
  
  char gnum[5]; //padded gradient direction id
  
  std::string curKey;
  std::string curVal;

  //take care of the baseline values
  for (unsigned int i = 0; i < num_baselines; i++)
    {
    sprintf(gnum, "%04d", i);
    curKey = "DWMRI_gradient_" + std::string(gnum);
    curVal = std::string("0 0 0");
    itk::EncapsulateMetaData<std::string>(newDictionary, curKey, curVal);
    }
  
  //take care of the rest of the gradient directions
  for (unsigned int i = 0; i < num_newgrads; i++)
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

template <class RealType>
//fills up a vnl_matrix with the gradients in the provided text file
void parseGradientFile(const std::string &gfname, vnl_matrix<RealType> &newgradients)
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

void getFiles( const std::string& sCaseFile, std::vector<std::string>& dwiFiles, std::vector<std::string>& deformationFiles, bool bNoDeformations )
{
  // parses the case file (sCaseFile) and returns all the file names
  // for the dwis and for the corresponding deformation files

  dwiFiles.clear();
  deformationFiles.clear();

  std::cout << "Reading: " << sCaseFile << std::endl;

  std::ifstream in( sCaseFile.c_str() );
  
  if ( in.is_open() ) 
    {

    std::string line;
    
    while ( getline ( in, line ) ) 
      {
      std::string::size_type iI = line.find_first_not_of ( " \t\n\v" );
      
      if ( iI != std::string::npos && line[iI] != '#'  )
	{
	// found something that needs to be processed, read it
	std::istringstream ins;
	std::string sub;
	
	ins.str( line );
	
	ins >> sub;
	dwiFiles.push_back( sub );
	std::cout << "DWI = " << sub << "; ";
	
	if ( !bNoDeformations )
	  {
	  ins >> sub;
	  deformationFiles.push_back( sub );
	  std::cout << "DF = " << sub << std::endl;
	  }
	
	}
      }
    } 
  else
    {
    std::cout << "ERROR: Could not open '" << sCaseFile << "' as input file. ABORT." << std::endl;
    exit( EXIT_FAILURE );
    }
}

}
