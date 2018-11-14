//
//  System.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "System.hpp"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace arma;

System::System()
{
    M = 0;
    Jx = 1.;
    Jy = 0.5;
    Jz = 1.;
}

void System::Add(Site s)
{
    Sites.push_back(s);
}

void System::Simulate()
{
    if(round(log2(Sites.size())) != log2(Sites.size()))
    {
        cout << "System::Simulate -> #sites must be a power of two.\n";
        exit(1);
    }

    //Initialize blocks with sites
    for(int i=0; i<Sites.size(); i++)
    {
        Blocks.push_back(Block(Sites[i]));
    }

    // Iterates until one block is left
    while(Blocks.size() != 1)
    {
        cout << "The size of the system:" << Blocks.size() << endl;
        Iterate();
    }

    //cout << Blocks[0].GetDissipator(0) << endl;
}

void System::Iterate()
{
    // TODO: if the efficiency is compromised, 
    //we suggest to adopt the following solution:
    //swap + erase the useless blocks.
    std::vector<Block> MergedBlocks;
    for(int i=0; i<Blocks.size(); i+=2)
    {
        cout << "Fin qui tutto bene for the " << i << "time\n";
        MergedBlocks.push_back(Merge(Blocks[i], Blocks[i+1]));
    }
    if (MergedBlocks.size() != Blocks.size()/2)
    {
        cout << "System::Iterate -> MergedBlocks and Blocks dimensions do not match.\n";
        exit(1);
    }
    Blocks = MergedBlocks;
}

Block System::Merge(const Block& b1, const Block& b2)
{
    if(b1.C.size() != b2.C.size())
    {
        cout << "System::Merge -> There is an incongruence in the dimensions of the Cs.\n";
    }

    cx_mat evectDM12;
    Block MergedBlock;

    // Merging
    MergedBlock.H = kron(b1.I, b2.H) + kron(b1.H, b2.I) - Jx*kron(b1.Sx, b2.Sx) - Jy*kron(b1.Sy, b2.Sy) - Jz*kron(b1.Sz, b2.Sz);
    cout << "Hdim: ";
    PrintSize(MergedBlock.H);
    
    for(int k=0; k<pow(2, b1.C.size()); k++)
    {
        if(k<pow(2, b1.C.size()-1))
        {
            MergedBlock.C.push_back(kron(b1.C[k], b2.I));
        }
        else
        {
            MergedBlock.C.push_back(kron(b1.I, b2.C[k-pow(2, b1.C.size()-1)]));
        }
    }

    cout << "Cdim: ";
    PrintSize(MergedBlock.C[0]);
    MergedBlock.Sx = kron(b1.Sx, b2.I);
    cout << "SxDim: ";
    PrintSize(MergedBlock.Sx);
    MergedBlock.Sy = kron(b1.Sy, b2.I);
    cout << "SyDim: ";
    PrintSize(MergedBlock.Sy);
    MergedBlock.Sz = kron(b1.Sz, b2.I);
    cout << "SzDim: ";
    PrintSize(MergedBlock.Sz);
    MergedBlock.I = kron(b1.I, b2.I);
    cout << "IDim: ";
    PrintSize(MergedBlock.I);

    // Renormalization
    evectDM12 = kron(b1.GetBlockEvectDM(), b2.GetBlockEvectDM());
    cout << "evectDM12: ";
    PrintSize(evectDM12);

    int Nkeep = std::min<int>(M, evectDM12.n_cols);
    cx_mat O = evectDM12.tail_cols(Nkeep);
    cout << "ODim: ";
    PrintSize(O);

    MergedBlock.H = trans(O)*MergedBlock.H*O;
    MergedBlock.Sx = trans(O)*MergedBlock.Sx*O;
    MergedBlock.Sy = trans(O)*MergedBlock.Sy*O;
    MergedBlock.Sz = trans(O)*MergedBlock.Sz*O;
    MergedBlock.I = trans(O)*MergedBlock.I*O;
    
    vector<cx_mat> Dslice;
    for(int k=0; k<pow(2, b1.C.size()); k++)
    {
        MergedBlock.C[k] = trans(O)*MergedBlock.C[k]*O;
        Dslice.push_back(0.5*(-kron(MergedBlock.C[k].st()*conj(MergedBlock.C[k]), MergedBlock.I)
                                - kron(MergedBlock.I, trans(MergedBlock.C[k])*MergedBlock.C[k])
                                 + 2.*kron(conj(MergedBlock.C[k]),MergedBlock.C[k])));
    }
        
    MergedBlock.D = 0.*Dslice[0];
    for(int k=0; k<pow(2, b1.C.size()); k++)
    {
        MergedBlock.D += Dslice[k];
    }

    return MergedBlock;
}

void System::SetCornerSize(int m)
{
    M = m;
}

void System::PrintSize(arma::cx_mat& m)
{
    cout << m.n_rows << "x" << m.n_cols << endl;
}