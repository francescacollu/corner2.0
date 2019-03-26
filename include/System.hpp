//
//  System.hpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#ifndef System_hpp
#define System_hpp

#include "Site.hpp"
#include "Block.hpp"
#include <vector>

class System
{
    // size of the corner-space
    int M;
    std::vector<Site> Sites;
    std::vector<Block> Blocks;

    // coupling spin coefficients
    float Jx;
    float Jy;
    float Jz;

    // Each iteration merges blocks two by two
    void Iterate();

    // Merging two blocks corresponding to indices i, j
    Block Merge(const Block& b1, const Block& b2, bool IsLeft);

public:
    // The method adds one site
    void Add(Site);

    // The heart of the program
    void Simulate();

    // Set size of the corner-space
    void SetCornerSize(int m);

    // Set the coupling constants
    void SetCouplingConstants(float jx, float jy, float jz);

    System();

    // Calculate the expectation value of an operator
    void GetExpValue(const char *file_name);

    //Calculate the two-point correlation function between two sites
    void Get2PCorrelationFunction(const char* file_CorrFunc);

    //Calculate the convergence with respect to M
    double Convergence(int m);

    //Calculate spin current
    void GetSpinCurrent(const char* file_SpinCurr);

private:
    static void PrintSize(arma::cx_mat&);
};

#endif /* System_hpp */
