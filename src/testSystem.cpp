//
//  testSystem.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "System.hpp"
#include <iostream>

using namespace std;

int main()
{
    System sy;
    
    sy.Add(Site());
    sy.Add(Site(Site::DissipatorType::ZUp));
    sy.Add(Site(Site::DissipatorType::ZDown));
    sy.Add(Site());

    sy.SetCornerSize(10);

    sy.Simulate();
    return 0;
}
