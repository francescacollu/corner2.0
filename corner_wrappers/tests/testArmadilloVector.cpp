#include "armadilloVector.h"
#include <iostream>

using namespace corner;

int main()
{
  ArmadilloVector v;
  v << 4 << 3;
  
  std::cout << "Size: " << v.size() << std::endl;
  std::cout << v.get(0) << '\t' << v.get(1) << std::endl;
  
  return 0;
}