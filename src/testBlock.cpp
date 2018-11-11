//
//  testBlock.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include <iostream>
#include "Block.hpp"
#include "Site.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    Site s;
    Block b = Block(Site(Site::DissipatorType::ZUp));
    
    cout << b.GetDissipator(0) << endl;
    
    return 0;
}
