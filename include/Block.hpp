//
//  Block.hpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#ifndef Block_hpp
#define Block_hpp

#include "Site.hpp"
#include "HMatrix.hpp"
#include <armadillo>

class Block
{
    std::vector<arma::cx_mat> C; // Dissipator
    std::vector<arma::cx_mat> sigmaZ;
    arma::cx_mat I;
    arma::cx_mat Sx;
    arma::cx_mat Sy;
    arma::cx_mat Sz;
    arma::cx_mat H;
    arma::cx_mat dm;
    arma::cx_mat L;
    arma::cx_mat D;
    arma::cx_mat evectDM;

    Block(Site);
    Block();

public:

    // System is a friend class of class Block
    friend class System;

    // Returns the i-th dissipator
    const arma::cx_mat& GetDissipator(int) const;

    // Returns the dm of the block
    const arma::cx_mat& GetBlockEvectDM() const;    
};

#endif /* Block_hpp */
