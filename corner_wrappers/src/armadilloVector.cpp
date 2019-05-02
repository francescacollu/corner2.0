#include "armadilloVector.h"
#include "corner_exceptions.h"

#include <iostream>

using namespace corner;

ArmadilloVector::ArmadilloVector(const arma::cx_vec& _v)
{
  v = _v;
}

double ArmadilloVector::norm() const
{
  return arma::norm(v);
}

std::shared_ptr<vector> ArmadilloVector::operator+(const vector& vec) const
{
  arma::cx_vec arma_vec;

  if (v.size() != vec.size())
    throw IncompatibleVectorException(*this, vec);

  arma_vec.resize(v.size());

  // Fill with vec
  for (int i = 0, len = vec.size(); i < len; i++)
  {
    arma_vec(i) = vec.get(i);
  }

  ArmadilloVector* new_vec = new ArmadilloVector(arma_vec);

  // Use armadillo built-in sum
  new_vec->v = v + arma_vec;

  return std::shared_ptr<ArmadilloVector>(new_vec);
}
std::shared_ptr<vector> ArmadilloVector::operator-(const vector&) const
{
  return nullptr;
}
std::shared_ptr<vector> ArmadilloVector::operator/(const std::complex<double>&) const
{
  return nullptr;
}
std::shared_ptr<vector> ArmadilloVector::operator*(const std::complex<double>&) const
{
  return nullptr;
}
size_t ArmadilloVector::size() const
{
  return v.n_rows;
}

vector& ArmadilloVector::operator<<(std::complex<double> x)
{
  v.resize(v.size()+1);
  v(v.size()-1) = x;

  return *this;
}

std::complex<double> ArmadilloVector::get(size_t index) const
{
  if (index > v.size()-1 || index < 0)
    throw std::out_of_range("ArmadilloVector::get out of range");

  std::complex<double> x = v(index);
  return x;
}

void ArmadilloVector::set(size_t index, std::complex<double> x)
{
  arma::cx_vec vec;
  vec.col(index) = arma::cx_double(x.real(), x.imag());
  return;
}
