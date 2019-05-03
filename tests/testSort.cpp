#include "eigenV.hpp"
#include <armadillo>
#include <iostream>

void print(const std::vector<corner::eigenV>& v, const char* description)
{
    std::cout << description << ": \n";
    for (std::vector<corner::eigenV>::const_iterator itr =  v.begin(); itr != v.end(); itr++)
    {
        std::cout << itr->eigenvalue << '\t';
    }
    std::cout << '\n';
}

int main()
{
    std::vector<corner::eigenV> eV;

    corner::eigenV pair1;
    pair1.eigenvalue = arma::cx_double(0,0);

    corner::eigenV pair2;
    pair2.eigenvalue = arma::cx_double(1, -2);
        
    corner::eigenV pair3;
    pair3.eigenvalue = arma::cx_double(-1, 2);

    corner::eigenV pair4;
    pair4.eigenvalue = arma::cx_double(1, 1);

    corner::eigenV pair5;
    pair5.eigenvalue = arma::cx_double(0, 1);

    corner::eigenV pair6;
    pair6.eigenvalue = arma::cx_double(2, 0);

    corner::eigenV pair7;
    pair7.eigenvalue = arma::cx_double(1, 2);


    eV.push_back(pair1);
    eV.push_back(pair2);
    eV.push_back(pair3);
    eV.push_back(pair5);
    eV.push_back(pair3);
    eV.push_back(pair1);
    eV.push_back(pair2);
    eV.push_back(pair7);
    eV.push_back(pair6);
    eV.push_back(pair4);
    eV.push_back(pair2);


    print(eV, "UNSORTED");

    std::sort(eV.begin(), eV.end(), corner::sortByAbs());

    print(eV, "PRE-GROUP");

    corner::groupWithSameAbs(eV);

    print(eV, "SORTED");
   
    return 0;
}
