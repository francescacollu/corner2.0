#include <armadillo>
#include <iostream>
#include <vector>

using namespace arma;

int main()
{
  cx_mat matrix = cx_mat().zeros(128*128, 128*128);
  std::cout << matrix.n_rows << "x" << matrix.n_cols << std::endl;

  std::vector<cx_mat> v;
  for (int i = 0; i < 3; i++)
    v.push_back(matrix);

  std::cout << "Ciaone: " << v.size() << std::endl;

  return 0;
}
