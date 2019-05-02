#pragma once

#include <armadillo>

namespace corner
{
    class eigenV
    {
        public:
        arma::cx_double eigenvalue;
        arma::cx_vec eigenvector;

        bool operator==(const eigenV& eV) const;
    };

    struct sortByReal
    {
        inline bool operator()(const eigenV& e1, const eigenV& e2) const
        {
            return (std::real(e1.eigenvalue) < std::real(e2.eigenvalue));
        }
    };

    struct sortByImag
    {
        inline bool operator()(const eigenV& e1, const eigenV& e2) const
        {
            return (std::imag(e1.eigenvalue) < std::imag(e2.eigenvalue));
        }
    };

    struct sortByAbs
    {
        inline bool operator()(const eigenV& e1, const eigenV& e2) const
        {
            return (std::abs(e1.eigenvalue) < std::abs(e2.eigenvalue));
        }
    };

    void groupWithSameAbs(std::vector<eigenV>& eV);
    void groupWithSameReal(std::vector<eigenV>::iterator beg, std::vector<eigenV>::iterator end);
    
}