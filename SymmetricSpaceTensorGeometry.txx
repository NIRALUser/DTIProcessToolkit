// -*- Mode: C++ -*-

template<class T, unsigned int dimension>
T 
SymmetricSpaceTensorGeometry<T, dimension>
::InnerProduct(const TensorType & base, const TangentType & v,
               const TangentType & w)
{
  EigenValuesArrayType eigenValues;
  EigenVectorsMatrixType eigenVectors;
  MatrixType diag, diagInv;
  vnl_matrix<T> g, gInv;
  vnl_matrix<T> gTrans, gInvTrans;
  MatrixType vMatrix, wMatrix;
  int i;

  OrientedEigensystem(base, eigenValues, eigenVectors);

  for(i = 0; i < dimension; i++)
  {
    if(eigenValues[i] <= 0.0)
      return 0; // Need to add error handling here!!

    diag[i][i] = sqrt(eigenValues[i]);
    diagInv[i][i] = (1.0 / diag[i][i]);
  }

  g = eigenVectors.GetTranspose() * diag.GetVnlMatrix();
  gTrans = g.transpose();

  gInv = diagInv.GetVnlMatrix() * eigenVectors.GetVnlMatrix();
  gInvTrans = gInv.transpose();

  TensorToMatrix(v, vMatrix);
  TensorToMatrix(w, wMatrix);

  return vnl_trace(gInv * vMatrix.GetVnlMatrix() * gInvTrans *
                   gInv * wMatrix.GetVnlMatrix() * gInvTrans);
}

template<class T, unsigned int dimension>
typename SymmetricSpaceTensorGeometry<T, dimension>::TensorType
SymmetricSpaceTensorGeometry<T, dimension>
::ExpMap(const TensorType & base, const TangentType & v)
{
  EigenValuesArrayType eigenValues;
  EigenVectorsMatrixType eigenVectors;
  MatrixType diag;
  MatrixType diagInv;
  vnl_matrix<T> g;
  vnl_matrix<T> gInv;
  MatrixType y;
  MatrixType result;
  MatrixType vMatrix;
  TensorType tensor;
  int i;

  OrientedEigensystem(base, eigenValues, eigenVectors);
  for(i = 0; i < dimension; i++)
  {
    if(eigenValues[i] <= 0.0)
      return tensor; // Need to add error handling here!!

    diag[i][i] = sqrt(eigenValues[i]);
    diagInv[i][i] = 1.0 / diag[i][i];
  }

  g = eigenVectors.GetTranspose() * diag.GetVnlMatrix();
  gInv = diagInv.GetVnlMatrix() * eigenVectors.GetVnlMatrix();

  TensorToMatrix(v, vMatrix);

  y = gInv * vMatrix.GetVnlMatrix() * gInv.transpose();

  MatrixToTensor(y, tensor);
  OrientedEigensystem(tensor, eigenValues, eigenVectors);

  for(i = 0; i < dimension; i++)
    diag[i][i] = exp(eigenValues[i]);

  result = (g * eigenVectors.GetTranspose() * diag.GetVnlMatrix() * 
            eigenVectors.GetVnlMatrix() * g.transpose());

  MatrixToTensor(result, tensor);
  return tensor;
}

template<class T, unsigned int dimension>
typename SymmetricSpaceTensorGeometry<T, dimension>::TangentType
SymmetricSpaceTensorGeometry<T, dimension>
::LogMap(const TensorType & base, const TensorType & p)
{
  EigenValuesArrayType eigenValues;
  EigenVectorsMatrixType eigenVectors;
  itk::Matrix<T, dimension> diag;
  itk::Matrix<T, dimension> diagInv;
  vnl_matrix<T> g;
  vnl_matrix<T> gInv;
  itk::Matrix<T, dimension> y;
  itk::Matrix<T, dimension> result;
  itk::Matrix<T, dimension> pMatrix;
  TensorType tensor;
  int i;

  OrientedEigensystem(base, eigenValues, eigenVectors);
  for(i = 0; i < dimension; i++)
  {
    if(eigenValues[i] <= 0.0)
      return tensor; // Need to add error handling here!!

    diag[i][i] = sqrt(eigenValues[i]);
    diagInv[i][i] = 1.0 / diag[i][i];
  }

  TensorToMatrix(p, pMatrix);

  g = eigenVectors.GetTranspose() * diag.GetVnlMatrix();
  gInv = diagInv.GetVnlMatrix() * eigenVectors.GetVnlMatrix();

  y = gInv * pMatrix.GetVnlMatrix() * gInv.transpose();

  MatrixToTensor(y, tensor);
  OrientedEigensystem(tensor, eigenValues, eigenVectors);

  for(i = 0; i < dimension; i++)
    diag[i][i] = log(eigenValues[i]);

  result = (g * eigenVectors.GetTranspose() * diag.GetVnlMatrix() * 
            eigenVectors.GetVnlMatrix() * g.transpose());

  MatrixToTensor(result, tensor);
  return tensor;
}

template<class T, unsigned int dimension>
void
SymmetricSpaceTensorGeometry<T, dimension>
::OrientedEigensystem(const TensorType & p, EigenValuesArrayType & eigenValues,
                      EigenVectorsMatrixType & eigenVectors)
{
  T tempVal;
  int i;

  p.ComputeEigenAnalysis(eigenValues, eigenVectors);

  // If det is negative, swap first two eigenvectors/values
  if(dimension >= 2 && vnl_det(eigenVectors.GetVnlMatrix()) < 0)
  {
    tempVal = eigenValues[0];
    eigenValues[0] = eigenValues[1];
    eigenValues[1] = tempVal;

    for(i = 0; i < dimension; i++)
    {
      tempVal = eigenVectors[0][i];
      eigenVectors[0][i] = eigenVectors[1][i];
      eigenVectors[1][i] = tempVal;
    }
  }
}

template<class T, unsigned int dimension>
typename SymmetricSpaceTensorGeometry<T, dimension>::TensorType
SymmetricSpaceTensorGeometry<T, dimension>
::GroupAction(const TensorType & p, const MatrixType & g)
{
  MatrixType resultMatrix;
  TensorType resultTensor;
  vnl_matrix<T> gVnl;
  vnl_matrix<T> pVnl;

  TensorToMatrix(p, resultMatrix);

  gVnl = g.GetVnlMatrix();
  pVnl = resultMatrix.GetVnlMatrix();

  resultMatrix = (gVnl * pVnl * gVnl.transpose());

  MatrixToTensor(resultMatrix, resultTensor);
  return resultTensor;
}

template<class T, unsigned int dimension>
void
SymmetricSpaceTensorGeometry<T, dimension>
::TensorToMatrix(const TensorType & p, MatrixType & m)
{
  int i, j, k;

  k = 0;
  for(i = 0; i < dimension; i++)
  {
    for(j = i; j < dimension; j++, k++)
    {
      m(i, j) = p[k];
      m(j, i) = p[k];
    }
  }
}

template<class T, unsigned int dimension>
void
SymmetricSpaceTensorGeometry<T, dimension>
::MatrixToTensor(const MatrixType & m, TensorType & p)
{
  int i, j, k;

  k = 0;
  for(i = 0; i < dimension; i++)
  {
    for(j = i; j < dimension; j++, k++)
    {
      p[k] = 0.5 * (m(i, j) + m(j, i));
    }
  }
}
