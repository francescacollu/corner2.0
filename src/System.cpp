//
//  System.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "System.hpp"
#include <cmath>

using namespace std;

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

    cout << Blocks[0].GetDissipator(0) << endl;
}

void System::Iterate()
{
    // TODO: if the efficiency is compromised, 
    //we suggest to adopt the following solution:
    //swap + erase of the useless blocks.
    std::vector<Block> MergedBlocks;
    for(int i=0; i<Blocks.size(); i+=2)
    {
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
    // FIXME: dummy Merge which substitutes the couple with the first block (twin thing)
    return b1;
}

void System::SetCornerSize(int m)
{
    M = m;
}



