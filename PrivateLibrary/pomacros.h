#ifndef POMACROS_H
#define POMACROS_H

#include <iostream>
#include <string>

#define verboseMessage(X)                       \
  if( VERBOSE ) std::cout << (X) << std::endl

void requireequal(bool test, const std::string & s, bool allowanyway)
{
  if( !test )
    {
    std::cerr << "WARNING: Trying to use two images with inconsistent meta-data" << std::endl;
    std::cerr << s << " are not equal." << std::endl;
    if( !allowanyway )
      {
      std::cerr << "Rerun with '--force' to run program anyway" << std::endl;
      exit(EXIT_FAILURE);
      }
    else
      {
      std::cerr << "WARNING: Overriding image consistency checks" << std::endl;
      std::cerr
        <<
      "The program may crash, use up all your memory, get stuck in an infinite loop, or something even more biazzare."
        << std::endl;
      }
    }
}

#endif
