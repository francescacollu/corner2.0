//
//  Block.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "Block.hpp"

using namespace arma;

Block::Block(Site s)
{
    // Initial fill: one element
    C.push_back(arma::cx_mat());
    
    if (s.GetDissipator() == Site::DissipatorType::Empty)
    {
        C[0] << cx_double(1., 0.) << cx_double(0., 0.) << endr << cx_double(0., 0.) << cx_double(1., 0.) << endr;
    }
    else if (s.GetDissipator() == Site::DissipatorType::ZUp)
    {
        C[0] << cx_double(0., 0.) << cx_double(1., 0.) << endr << cx_double(0., 0.) << cx_double(0., 0.) << endr;
    }
    else if (s.GetDissipator() == Site::DissipatorType::ZDown)
    {
        C[0] << cx_double(0., 0.) << cx_double(0., 0.) << endr << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    }
    else
    {
        std::cout << "Block::Block -> Type has not been recognized.\n";
        exit(1);
    }
}

const arma::cx_mat& Block::GetDissipator(int index) const
{
    if (index > C.size() || index < 0)
    {
        cout << "Block::GetDissipator() -> index out of bounds.\n";
    }
    return C[index];
}

