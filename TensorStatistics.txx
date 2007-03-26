// -*- Mode: C++ -*-

template<class T, unsigned int dimension>
void
TensorStatistics<T, dimension>
::ComputeMean(const TensorListPointerType tensorList, TensorType & mean) const
{
  TangentType tangent(0.0), initialTangent;
  T lastNormSquared, normSquared, initialNormSquared;
  T currStepSize = stepSize;
  int i;

  if(tensorList->Size() == 0)
    return;

  mean = tensorList->ElementAt(0);

  for(i = 0; i < tensorList->Size(); i++)
    tangent += tensGeometry->LogMap(mean, tensorList->ElementAt(i));

  tangent = tangent * (1.0 / ((T) tensorList->Size()));
  initialTangent = tangent;

  initialNormSquared = tensGeometry->NormSquared(mean, tangent);
  lastNormSquared = initialNormSquared;
  while(lastNormSquared >= EPSILON)
  {
    mean = tensGeometry->ExpMap(mean, tangent * currStepSize);

    tangent.Fill(0.0);
    for(i = 0; i < tensorList->Size(); i++)
      tangent += tensGeometry->LogMap(mean, tensorList->ElementAt(i));

    tangent = tangent * (1.0 / ((T) tensorList->Size()));

    normSquared = tensGeometry->NormSquared(mean, tangent);
    if(normSquared >= lastNormSquared)
    {
      currStepSize *= 0.5;
      mean = tensorList->ElementAt(0);
      lastNormSquared = initialNormSquared;
      tangent = initialTangent;
    }
    else
      lastNormSquared = normSquared;
  }
}

template<class T, unsigned int dimension>
void
TensorStatistics<T, dimension>
::ComputeWeightedAve(const ScalarListPointerType weightList,
                     const TensorListPointerType tensorList,
                     TensorType & weightedAve) const
{
  TangentType tangent(0.0), initialTangent;
  T lastNormSquared, normSquared, initialNormSquared;
  T currStepSize = stepSize;
  int i;

  if(tensorList->Size() == 0 ||
     tensorList->Size() != weightList->Size())
    return;

  weightedAve = tensorList->ElementAt(0);

  for(i = 0; i < tensorList->Size(); i++)
    tangent += weightList->ElementAt(i) *
      tensGeometry->LogMap(weightedAve, tensorList->ElementAt(i));

  initialTangent = tangent;

  initialNormSquared = tensGeometry->NormSquared(weightedAve, tangent);
  lastNormSquared = initialNormSquared;
  while(lastNormSquared >= EPSILON)
  {
    weightedAve = tensGeometry->ExpMap(weightedAve, tangent * currStepSize);

    tangent.Fill(0.0);
    for(i = 0; i < tensorList->Size(); i++)
      tangent += weightList->ElementAt(i) *
        tensGeometry->LogMap(weightedAve, tensorList->ElementAt(i));

    normSquared = tensGeometry->NormSquared(weightedAve, tangent);
    if(normSquared >= lastNormSquared)
    {
      currStepSize *= 0.5;
      weightedAve = tensorList->ElementAt(0);
      lastNormSquared = initialNormSquared;
      tangent = initialTangent;
    }
    else
      lastNormSquared = normSquared;
  }
}

template<class T, unsigned int dimension>
void
TensorStatistics<T, dimension>
::ComputeMeanAndCovariance(const TensorListPointerType tensorList,
                           TensorType & mean,
                           CovarianceType & covariance) const
{
  const int size = dimension * (dimension + 1) / 2;
  TensorType logTens;
  int i, j, k, index;

  ComputeMean(tensorList, mean);
  covariance.Fill(0.0);

  for(i = 0; i < tensorList->Size(); i++)
  {
    logTens = tensGeometry->LogMap(mean, tensorList->ElementAt(i));

    index = 0;
    for(j = 0; j < size; j++)
    {
      for(k = j; k < size; k++, index++)
      {
        covariance[index] += logTens[j] * logTens[k];
      }
    }
  }

  covariance = covariance * (1.0 / ((T) tensorList->Size() + 1));
}

template<class T, unsigned int dimension>
void
TensorStatistics<T, dimension>
::ComputeMeanAndPGA(const TensorListPointerType tensorList, TensorType & mean,
                    PGAVariancesArrayType & pgaVariances,
                    PGAVectorsMatrixType & pgaVectorsMatrix) const
{
  CovarianceType covariance;

  ComputeMeanAndCovariance(tensorList, mean, covariance);

  covariance.ComputeEigenAnalysis(pgaVariances, pgaVectorsMatrix);
}

template<class T, unsigned int dimension>
typename TensorStatistics<T, dimension>::TensorType
TensorStatistics<T, dimension>
::RandomGaussianTensor(const TensorType & mean, T sigma) const
{
  double x, y, r2;
  int i, dim;
  TensorType logResult(0.0);

  // Generate random gaussians with polar Box-Muller method
  dim = (dimension * (dimension + 1)) / 2;
  for(i = 0; i < dim; i++)
  {
    r2 = 0;
    while(r2 > 1.0 || r2 == 0)
    {
      x = (2.0 * (double) rand()) / (double)(RAND_MAX + 1.0) - 1.0;
      y = (2.0 * (double) rand()) / (double)(RAND_MAX + 1.0) - 1.0;

      r2 = x * x + y * y;
    }

    logResult[i] = sigma * y * sqrt(-2.0 * log(r2) / r2);
  }

  return tensGeometry->ExpMap(mean, logResult);
}
