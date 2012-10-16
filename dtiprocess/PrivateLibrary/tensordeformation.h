#ifndef TENSORDEFORMATION_H
#define TENSORDEFORMATION_H

#include "dtitypes.h"

// warping functions
TensorImageType::Pointer createROT(TensorImageType::Pointer, const std::string &,int doffiletype);
TensorImageType::Pointer createWarp(TensorImageType::Pointer, 
                                    DeformationImageType::Pointer, 
                                    TensorReorientationType, InterpolationType);


#endif
