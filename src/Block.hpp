//
//  Block.hpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#ifndef Block_hpp
#define Block_hpp

#include "Site.hpp"
#include <armadillo>

class Block
{
    std::vector<arma::cx_mat> C; // Dissipator
public:
    Block(Site);
    // Returns the i-th dissipator
    const arma::cx_mat& GetDissipator(int) const;
    
};

#endif /* Block_hpp */
