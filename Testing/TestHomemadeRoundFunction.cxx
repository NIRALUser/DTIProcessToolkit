#include "itkDiffusionTensor3DReconstructionImageFilterBase.h"
#include <vector>

int main(int argc, char* argv[])
{
	/*
	Negative numbers are not tested because the homemade "round()"
	function is not implemented for negative numbers because it
	is only used to round numbers from the exponential function:
	See Library/itkDiffusionTensor3DReconstructionImageFilterBase.txx:156:
	bit.Set(static_cast<GradientPixelType>(round(exp(D[6]))));
	*/

	std::vector<double> vec;
	std::vector<long> vecResult;
	vec.push_back(0.1);
	vecResult.push_back(0);

	vec.push_back(0.2);
	vecResult.push_back(0);

	vec.push_back(0.5);
	vecResult.push_back(1);

	vec.push_back(1.1);
	vecResult.push_back(1);

	vec.push_back(2.71828);
	vecResult.push_back(3);

	vec.push_back(15.236);
	vecResult.push_back(15);

	vec.push_back(5488);
	vecResult.push_back(5488);

	vec.push_back(253.0);
	vecResult.push_back(253);

	vec.push_back(154.45);
	vecResult.push_back(154);

	for(unsigned int i=0;i<vec.size();i++)
	{
		long rounded = round(vec[i]);
		if(rounded!=vecResult[i]) return -1;
//		std::cout<<"round("<<vec[i]<<") = "<<rounded<<std::endl;
	}

	return 0;
}
