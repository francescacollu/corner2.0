#pragma once
#include "vector.h"
#include <armadillo>
#include <complex>

namespace corner
{

  class ArmadilloVector:public vector
  {
  public:
    ArmadilloVector(const arma::cx_vec& _v);
    std::shared_ptr<vector> operator+(const vector&) const override;
    std::shared_ptr<vector> operator-(const vector&) const override;
    std::shared_ptr<vector> cdot(const vector&) const override;
    std::shared_ptr<vector> operator/(const std::complex<double>&) const override;
    std::shared_ptr<vector> operator*(const std::complex<double>&) const override;
    double norm() const override;
    size_t size() const override;
    vector& operator<<(std::complex<double>) override;
    std::complex<double> get(size_t index) const override;
    void set(size_t index, std::complex<double>) override;
    
  private:
    arma::cx_vec v;
  };
  
}
