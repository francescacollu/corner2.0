#pragma once


namespace corner
{

  class matrix
  {
  public:
    virtual eigenV GetEigen() = 0;
    virtual matrix* operator*(const matrix&) const = 0;
    virtual 
  };
  
}