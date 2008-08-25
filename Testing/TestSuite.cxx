#include <iostream>
#include "localTestMain.h"

int CommandLineTest(int argc, char* argv[])
{
  if(argc == 0)
    std::cerr << "No command line test specified" << std::endl;
      
  std::string command;
  for(int i = 1 ;  i <argc; ++i)
    command += std::string(argv[i]) + " ";

  std::cout << "Executing: " << std::endl;
  std::cout << command << std::endl;
  return system(command.c_str());
}

void RegisterTests()
{
  REGISTER_TEST(BesselTest );
  REGISTER_TEST(ScalarDeformationApplyTest ); 
  REGISTER_TEST(WarpFiberTest );
  REGISTER_TEST(WarpTensorTest );
  REGISTER_TEST(FiberIOTest );
  REGISTER_TEST(LLSTensorEstimateTest );
  REGISTER_TEST(WLSTensorEstimateTest );
  REGISTER_TEST(NLSTensorEstimateTest );
  REGISTER_TEST(CommandLineTest );
}
