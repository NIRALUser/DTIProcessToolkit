#include "itkDiffusionTensor3DReconstructionImageFilterBase.h"

int main(int , char* [])
{
	/*
	Negative numbers are not tested because the homemade "round()"
	function is not implemented for negative numbers because it
	is only used to round numbers from the exponential function:
	See Library/itkDiffusionTensor3DReconstructionImageFilterBase.txx:156:
	bit.Set(static_cast<GradientPixelType>(round(exp(D[6]))));
	*/
	double vec [] = { 0.1 , 0.2 , 0.5 , 1.1 , 2.71828 , 15.236 , 5488.0 , 253.0 , 154.45 } ;
	long vecResult[] = { 0 , 0 , 1 , 1 , 3 , 15 , 5488 , 253 , 154 } ;
	for( unsigned int i = 0 ; i < 9 ; i++ )
	{
		long rounded = round( vec[ i ] ) ;
		if( rounded != vecResult[ i ] )
		{
			std::cout << rounded << " " << vecResult[i] << std::endl ;
			return -1;
		}
//		std::cout<<"round("<<vec[i]<<") = "<<rounded<<std::endl;
	}

	return 0;
}
