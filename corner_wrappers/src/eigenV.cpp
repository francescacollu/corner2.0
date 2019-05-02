#include <eigenV.h>

using namespace corner;

bool eigenV::operator<(eigenV eV)
{
  if(abs(eigenvalue) < abs(eV.eigenvalue))
    return true;
  
  return false;
}

bool eigenV::operator>(eigenV eV)
{
  if(!(*this < eV) && !(*this == eV))
    return true;
  
  return false;
}

bool eigenV::operator==(eigenV eV)
{
  if(abs(eigenvalue) == abs(eV.eigenvalue) && eigenvector == eV.eigenvector)
    return true;
  
  return false;
}

bool eigenV::operator!=(eigenV eV)
{
  return !(*this == eV);
}