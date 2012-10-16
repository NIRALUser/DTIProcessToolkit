#ifndef _AUX_FUNCTIONS_
#define _AUX_FUNCTIONS_

#include <itkMetaDataObject.h>
#include <vector>
#include <string>
#include <itkArray2D.h>

namespace aux
{

template <class RealType>
void ConstructOutputMetaDataDictionary( itk::MetaDataDictionary & newDictionary, const itk::MetaDataDictionary & dictionary, unsigned int num_baselines, const vnl_matrix<RealType> &newgrads  );

template <class RealType>
void parseGradientFile(const std::string &gfname, vnl_matrix<RealType> &newgradients);

int getNumOutputGradients(const std::string &gfname);
void getFiles( const std::string& sCaseFile, std::vector<std::string>& dwiFiles, std::vector<std::string>& deformationFiles, bool bNoDeformations );

}

#include "auxFunctions.txx"

#endif
