#ifndef TENSORDEFORMATION_H
#define TENSORDEFORMATION_H

#include "dtitypes.h"

// warping functions
TensorImageType::Pointer createROT(TensorImageType::Pointer, const std::string &);
TensorImageType::Pointer createWarp(TensorImageType::Pointer, 
                                    DeformationImageType::Pointer, 
                                    TensorReorientationType, InterpolationType);


#endif
