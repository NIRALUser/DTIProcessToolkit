#ifndef TENSORDEFORMATION_H
#define TENSORDEFORMATION_H

#include "dtitypes.h"

// warping functions
void createROT(TensorImageType::Pointer, const std::string &, const std::string &);
void createWarp(TensorImageType::Pointer, const std::string &, const std::string &, const std::string &, TensorReorientationType, InterpolationType);


#endif
