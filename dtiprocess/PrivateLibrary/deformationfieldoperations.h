#ifndef DEFORMATIONFIELDIO_H
#define DEFORMATIONFIELDIO_H

#include "dtitypes.h"
#include <string>

enum DeformationFieldType {HField, Displacement};

DeformationImageType::Pointer readDeformationField(std::string warpfile, DeformationFieldType dft);


#endif
