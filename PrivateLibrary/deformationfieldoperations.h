#ifndef DEFORMATIONFIELDOPERATIONS_H
#define DEFORMATIONFIELDOPERATIONS_H

#include "dtitypes.h"
#include <string>

enum DeformationFieldType {HField, Displacement};

DeformationImageType::Pointer readDeformationField(std::string warpfile, DeformationFieldType dft);


#endif
