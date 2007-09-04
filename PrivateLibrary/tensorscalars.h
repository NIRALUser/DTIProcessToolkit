#ifndef TENSORSCALARS_H
#define TENSORSCALARS_H

#include "dtitypes.h"

// derived output functions
void createFA(TensorImageType::Pointer, const std::string &, bool);
void createFAGradient(TensorImageType::Pointer, const std::string &,double);
void createFAGradMag(TensorImageType::Pointer, const std::string &,double);
void createColorFA(TensorImageType::Pointer, const std::string &);
void createPrincipalEigenvector(TensorImageType::Pointer, const std::string &, const std::string &, GradientListType::Pointer);
void createMD(TensorImageType::Pointer, const std::string &, bool);
void createNegativeEigenValueLabel(TensorImageType::Pointer, const std::string &);

#endif
