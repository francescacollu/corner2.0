#pragma once

#include <memory>
#include <complex>

namespace corner
{

  class vector
  {
  public:
    virtual std::shared_ptr<vector> operator+(const vector&) const = 0;
    virtual std::shared_ptr<vector> operator-(const vector&) const = 0;
    virtual std::shared_ptr<vector> cdot(const vector&) const = 0;
    virtual std::shared_ptr<vector> operator/(const std::complex<double>&) const = 0;
    virtual std::shared_ptr<vector> operator*(const std::complex<double>&) const = 0;
    virtual double norm() const = 0;
    virtual size_t size() const = 0;
    virtual vector& operator<<(std::complex<double>) = 0;
    virtual std::complex<double> get(size_t) const = 0;
    virtual void set(size_t, std::complex<double>) = 0;
  };
  
}