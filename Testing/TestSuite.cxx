#include <iostream>
#include "localTestMain.h"

int CommandLineTest(int argc, char* argv[])
{
  if(argc == 0)
    std::cerr << "No command line test specified" << std::endl;
      
  std::string command;
  for(int i = 1 ;  i <argc; ++i)
    command += std::string(argv[i]) + " ";

  return system(command.c_str());
}

void RegisterTests()
{
  REGISTER_TEST(BesselTest );
  REGISTER_TEST(ScalarDeformationApplyTest ); 
//  REGISTER_TEST(WarpFiberTest );
//  REGISTER_TEST(WarpTensorTest );
  REGISTER_TEST(FiberIOTest );
  REGISTER_TEST(TensorEstimateTest );
  REGISTER_TEST(CommandLineTest );
}
