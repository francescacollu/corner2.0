#pragma once

#include <complex>
#include "vector.h"

namespace corner
{

  class eigenV
  {
  public:
    complex<double> eigenvalue;
    std::vector<corner::vector> eigenvectors;

    bool operator<(eigenV eV);
    bool operator>(eigenV eV);
    bool operator==(eigenV eV);
    bool operator!=(eigenV eV);
  };

}