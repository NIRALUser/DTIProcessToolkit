#ifndef TENSORSCALARS_H
#define TENSORSCALARS_H

#include "dtitypes.h"
#include <itkImage.h>

// derived output functions
template<class T>
typename itk::Image<T, 3>::Pointer createFA(TensorImageType::Pointer);

template<class T>
typename itk::Image<T, 3>::Pointer createMD(TensorImageType::Pointer);

template<class T>
typename itk::Image<T, 3>::Pointer createFro(TensorImageType::Pointer);

GradientImageType::Pointer createFAGradient(TensorImageType::Pointer, double);
RealImageType::Pointer createFAGradMag(TensorImageType::Pointer, double);
RGBImageType::Pointer createColorFA(TensorImageType::Pointer);
LabelImageType::Pointer createNegativeEigenValueLabel(TensorImageType::Pointer);

// TODO: remove filename
GradientImageType::Pointer createPrincipalEigenvector(TensorImageType::Pointer);

#endif
