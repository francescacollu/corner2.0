//
//  System.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "System.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>

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
        cout << "\n\nThe system contains " << Blocks.size() << " blocks.\n";
        Iterate();
    }
}

void System::Iterate()
{
    // TODO: if the efficiency is compromised, 
    //we suggest to adopt the following solution:
    //swap + erase the useless blocks.
    std::vector<Block> MergedBlocks;
    for(int i=0; i<Blocks.size(); i+=2)
    {
        cout << "Dealing with the merge " << i << " and " << i+1 << " of the chain of " << Blocks.size() << " blocks.\n";
        if(i%4 == 0)
        {
            MergedBlocks.push_back(Merge(Blocks[i], Blocks[i+1], true));
        }
        else 
        {
            MergedBlocks.push_back(Merge(Blocks[i], Blocks[i+1], false));
        }
    }
    if (MergedBlocks.size() != Blocks.size()/2)
    {
        cout << "System::Iterate -> MergedBlocks and Blocks dimensions do not match.\n";
        exit(1);
    }
    Blocks = MergedBlocks;
}

Block System::Merge(const Block& b1, const Block& b2, bool IsLeft)
{
    if(b1.C.size() != b2.C.size())
    {
        cout << "System::Merge -> There is an incongruence in the dimensions of the Cs.\n";
    }

    cx_mat evectDM12;
    Block MergedBlock;

    // Merging
    MergedBlock.H = kron(b1.I, b2.H) + kron(b1.H, b2.I) - Jx*kron(b1.Sx, b2.Sx) - Jy*kron(b1.Sy, b2.Sy) - Jz*kron(b1.Sz, b2.Sz);
    
    for(int k=0; k<2*b1.C.size(); k++)
    {
        if(k<b1.C.size())
        {
            MergedBlock.C.push_back(kron(b1.C[k], b2.I));
            MergedBlock.sigmaZ.push_back(kron(b1.sigmaZ[k], b2.I));
        }
        else
        {
            MergedBlock.C.push_back(kron(b1.I, b2.C[k-b1.C.size()]));
            MergedBlock.sigmaZ.push_back(kron(b1.I, b2.sigmaZ[k-b1.sigmaZ.size()]));
        }
    }

    if(IsLeft)
    {
        MergedBlock.Sx = kron(b1.I, b2.Sx);
        MergedBlock.Sy = kron(b1.I, b2.Sy);
        MergedBlock.Sz = kron(b1.I, b2.Sz);
        MergedBlock.I = kron(b1.I, b2.I);
    }
    else
    {
        MergedBlock.Sx = kron(b1.Sx, b2.I);
        MergedBlock.Sy = kron(b1.Sy, b2.I);
        MergedBlock.Sz = kron(b1.Sz, b2.I);
        MergedBlock.I = kron(b1.I, b2.I);
    }

    // Renormalization
    cx_mat dm12 = kron(b1.dm, b2.dm);
    HMatrix DM12(dm12);
    evectDM12 = DM12.GetONBasis();

    int Nkeep = std::min<int>(M, evectDM12.n_cols);
    cx_mat O = evectDM12.tail_cols(Nkeep);

    MergedBlock.H = trans(O)*MergedBlock.H*O;
    MergedBlock.Sx = trans(O)*MergedBlock.Sx*O;
    MergedBlock.Sy = trans(O)*MergedBlock.Sy*O;
    MergedBlock.Sz = trans(O)*MergedBlock.Sz*O;
    MergedBlock.I = trans(O)*MergedBlock.I*O;
    
    vector<cx_mat> Dslice;
    for(int k=0; k<2*b1.C.size(); k++)
    {
        MergedBlock.sigmaZ[k] = trans(O)*MergedBlock.sigmaZ[k]*O;
        MergedBlock.C[k] = trans(O)*MergedBlock.C[k]*O;
        Dslice.push_back(0.5*(-kron(MergedBlock.C[k].st()*conj(MergedBlock.C[k]), MergedBlock.I)
                                - kron(MergedBlock.I, trans(MergedBlock.C[k])*MergedBlock.C[k])
                                + 2.*kron(conj(MergedBlock.C[k]),MergedBlock.C[k])));
    }

    MergedBlock.D = 0.*Dslice[0];
    for(int k=0; k<2*b1.C.size(); k++)
    {
        MergedBlock.D += Dslice[k];
    }

    MergedBlock.L = cx_double(0., 1.)*kron(MergedBlock.H.st(), MergedBlock.I) - cx_double(0., 1.)*kron(MergedBlock.I, MergedBlock.H) + MergedBlock.D;

    HMatrix liouv(MergedBlock.L);
    MergedBlock.dm = liouv.GetSteadyStateDM();
    HMatrix DM(MergedBlock.dm);
    MergedBlock.evectDM = DM.GetONBasis();

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

void System::GetExpValue()
{
    if(Blocks.size() != 1)
    {
        cout << "System::GetExpValue -> The merge hasn't worked.";
    }

    ofstream myfile("output.txt");

    for(int i=0; i<Sites.size(); i++)
    {
        myfile << i << "\t" << real(arma::trace(Blocks[0].sigmaZ[i]*Blocks[0].dm)) << endl;
        cout << "<sigmaZ[" << i << "]> = " << arma::trace(Blocks[0].sigmaZ[i]*Blocks[0].dm) << endl;
    }
}