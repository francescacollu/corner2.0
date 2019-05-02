#include <armadillo>
#include <iostream>

using namespace arma;

int main()
{
    // Only test this if arpack should be available
    #ifdef USE_ARPACK
    arma::sp_cx_mat S;

    S.eye(65536, 65536);
    cx_vec eigval;
    cx_mat eigvec;

    eigs_gen(eigval, eigvec, S, S.n_rows);

    std::cout << eigval << std::endl;
    
    #endif

    return 0;
}
