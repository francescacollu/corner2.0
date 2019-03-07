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
#include "utility.hpp"

using namespace std;
using namespace arma;
using namespace corner;

System::System()
{
    M = 0;
    Jx = 1.;
    Jy = 1.;
    Jz = 1.;
}

void System::Add(Site s)
{
    Sites.push_back(s);
}

void System::Simulate()
{
    cout << __LINE__ << endl;
    if(round(log2(Sites.size())) != log2(Sites.size()))
    {
        cout << "System::Simulate -> #sites must be a power of two.\n";
        exit(1);
    }
    cout << __LINE__ << endl;
    //Initialize blocks with sites
    for(int i=0; i<Sites.size(); i++)
    {
        Blocks.push_back(Block(Sites[i]));
    }

    cout << __LINE__ << endl;
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
    MergedBlock.H = kron(b1.I, b2.H) + kron(b1.H, b2.I) + Jx*kron(b1.Sx, b2.Sx) + Jy*kron(b1.Sy, b2.Sy) + Jz*kron(b1.Sz, b2.Sz);
    
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
    cout << __LINE__ << endl;
    if(!DM12.IsHermitian()){cout << "Not hermitian\n";}
    if(!DM12.TraceOne()){cout << "Not trace 1 = " << trace(dm12) <<"\n";}
    if(!DM12.EigSumIsOne()){cout << "Sum not one\n";}
    if(!DM12.EigRealPositive()){cout << "eigenvalues not real positive:" << endl;}
    //check(DM12.IsDM(), "System::Merge", "This is not a DM");
    
    evectDM12 = DM12.GetONBasis();
    cout << __LINE__ << endl;

    int Nkeep = std::min<int>(M, evectDM12.n_cols);
    cx_mat O = evectDM12.tail_cols(Nkeep);

    MergedBlock.H = trans(O)*MergedBlock.H*O;
    MergedBlock.Sx = trans(O)*MergedBlock.Sx*O;
    MergedBlock.Sy = trans(O)*MergedBlock.Sy*O;
    MergedBlock.Sz = trans(O)*MergedBlock.Sz*O;
    MergedBlock.I = trans(O)*MergedBlock.I*O;
    
    cout << __LINE__ << endl;
    vector<cx_mat> Dslice;
    for(int k=0; k<2*b1.C.size(); k++)
    {
        MergedBlock.sigmaZ[k] = trans(O)*MergedBlock.sigmaZ[k]*O;
        MergedBlock.C[k] = trans(O)*MergedBlock.C[k]*O;
        Dslice.push_back(0.5*(-kron(MergedBlock.C[k].st()*conj(MergedBlock.C[k]), MergedBlock.I)
                                - kron(MergedBlock.I, trans(MergedBlock.C[k])*MergedBlock.C[k])
                                + 2.*kron(conj(MergedBlock.C[k]),MergedBlock.C[k])));
    }

    cout << __LINE__ << endl;
    MergedBlock.D = 0.*Dslice[0];
    for(int k=0; k<2*b1.C.size(); k++)
    {
        MergedBlock.D += Dslice[k];
    }

    MergedBlock.L = cx_double(0., 1.)*kron(MergedBlock.H.st(), MergedBlock.I) - cx_double(0., 1.)*kron(MergedBlock.I, MergedBlock.H) + MergedBlock.D;

    cout << __LINE__ << endl;
    HMatrix liouv(MergedBlock.L, true);
    cout << __LINE__ << endl;
    MergedBlock.dm = liouv.GetSteadyStateDM();
    cout << __LINE__ << endl;
    HMatrix DM(MergedBlock.dm);
    cout << __LINE__ << endl;
    MergedBlock.evectDM = DM.GetONBasis();

    return MergedBlock;
}

void System::SetCornerSize(int m)
{
    M = m;
}

void System::SetCouplingConstants(float jx, float jy, float jz)
{
    Jx = jx;
    Jy = jy;
    Jz = jz;
}

void System::PrintSize(arma::cx_mat& m)
{
    cout << m.n_rows << "x" << m.n_cols << endl;
}

void System::GetExpValue(const char *file_name)
{
    if(Blocks.size() != 1)
    {
        cout << "System::GetExpValue -> The merge hasn't worked.";
    }

    ofstream myfile(file_name);

    for(int i=0; i<Sites.size(); i++)
    {
        myfile << i << "\t" << real(arma::trace(Blocks[0].sigmaZ[i]*Blocks[0].dm)) << endl;
        cout << "<sigmaZ[" << i << "]> = " << arma::trace(Blocks[0].sigmaZ[i]*Blocks[0].dm) << endl;
    }
}

void System::Get2PCorrelationFunction(const char* file_CorrFunc, int i)
{
    ofstream myfile(file_CorrFunc);
    double tr0 = real(arma::trace(Blocks[0].sigmaZ[i]*Blocks[0].dm));
    
    for(int j=0; j<Sites.size(); j++)
    {
        double tr1 = real(arma::trace(Blocks[0].sigmaZ[j]*Blocks[0].dm));
        double tr01 = real(arma::trace(Blocks[0].sigmaZ[i]*Blocks[0].sigmaZ[j]*Blocks[0].dm));
        myfile << j << "\t" << tr01-tr0*tr1 << endl;
        cout << "tr01 = " << tr01 << "\t tr0= " << tr0 << "\t tr1 = " << tr1 << endl;
    }

}

double System::Convergence(int m)
{
    cout << m << "\t" << real(arma::trace(Blocks[0].sigmaZ[0]*Blocks[0].dm)) << endl;  
    return real(arma::trace(Blocks[0].sigmaZ[0]*Blocks[0].dm));  
}

// cx_double System::SpinCurrent(const char* file_SpinCurr, int i)
// {
//    ofstream myfile(file_SpinCurr);
//     for(int j=0; j<Sites.size(); j++)
//     {
//         double sc = real(arma::trace(Blocks[0].sigmaZ[j]**Blocks[0].dm));
//         myfile << j << "\t" << sc << endl;
//         cout << j << "\t" << sc << endl;
//     }

// }