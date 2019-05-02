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
    // Initial fill of dissipator: one element
    C.push_back(arma::cx_mat());
    sigmaZ.push_back(arma::cx_mat());
    sigmaX.push_back(arma::cx_mat());
    sigmaY.push_back(arma::cx_mat());

    I << cx_double(1., 0.) << cx_double(0., 0.) << endr 
      << cx_double(0., 0.) << cx_double(1., 0.) << endr;

    Sx << cx_double(0., 0.) << cx_double(1., 0.) << endr 
       << cx_double(1., 0.) << cx_double(0., 0.) << endr;

    Sy << cx_double(0., 0.) << cx_double(0., -1.) << endr 
       << cx_double(0., 1.) << cx_double(0., 0.) << endr;

    Sz << cx_double(1., 0.) << cx_double(0., 0.) << endr 
       << cx_double(0., 0.) << cx_double(-1., 0.) << endr;

    H << cx_double(0., 0.) << cx_double(0., 0.) << endr 
      << cx_double(0., 0.) << cx_double(0., 0.) << endr;

    sigmaZ[0] << cx_double(1., 0.) << cx_double(0., 0.) << endr 
              << cx_double(0., 0.) << cx_double(-1., 0.) << endr;
    
    sigmaX[0] << cx_double(0., 0.) << cx_double(1., 0.) << endr 
              << cx_double(1., 0.) << cx_double(0., 0.) << endr;

    sigmaY[0] << cx_double(0., 0.) << cx_double(0., -1.) << endr 
              << cx_double(0., 1.) << cx_double(0., 0.) << endr;

    
    if (s.GetDissipator() == Empty)
    {
        C[0] << cx_double(1., 0.) << cx_double(0., 0.) << endr 
             << cx_double(0., 0.) << cx_double(1., 0.) << endr;
    }
    else if (s.GetDissipator() == ZUp)
    {
        C[0] << cx_double(0., 0.) << cx_double(1., 0.) << endr 
             << cx_double(0., 0.) << cx_double(0., 0.) << endr;
    }
    else if (s.GetDissipator() == epsZUp)
    {
        C[0] << cx_double(0., 0.) << cx_double(0.01, 0.) << endr 
             << cx_double(0., 0.) << cx_double(0., 0.) << endr;
    }
    else if (s.GetDissipator() == ZDown)
    {
        C[0] << cx_double(0., 0.) << cx_double(0., 0.) << endr 
             << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    }
    else if (s.GetDissipator() == epsZDown)
    {
        C[0] << cx_double(0., 0.) << cx_double(0., 0.) << endr 
             << cx_double(0.01, 0.) << cx_double(0., 0.) << endr;
    }
    else
    {
        std::cout << "Block::Block -> Dissipator type has not been recognized.\n";
        exit(1);
    }

    D = -0.5*kron(C[0].st()*conj(C[0]), I) - 0.5*kron(I, trans(C[0])*C[0]) + kron(conj(C[0]), C[0]);
    L = cx_double(0., 1.)*kron(H.st(), I) - cx_double(0., 1.)*kron(I, H) + D;
    HMatrix liouv(L, true);
    //cout << "eigenvalues of L: " << liouv.GetEigenvalues() << endl;
    dm = liouv.GetSteadyStateDM();
    HMatrix DM(dm);
    evectDM = DM.GetONBasis();
}

Block::Block()
{

}

const arma::cx_mat& Block::GetDissipator(int index) const
{
    if (index > C.size() || index < 0)
    {
        cout << "Block::GetDissipator() -> index out of bounds.\n";
    }
    return C[index];
}

const arma::cx_mat& Block::GetBlockEvectDM() const
{
    return evectDM;
}


