//
//  testSite.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include <iostream>
#include "Site.hpp"

using namespace std;

int main()
{
    Site sEmpty;
    Site sUp(Site::DissipatorType::ZUp);
    Site sDown(Site::DissipatorType::ZDown);
    
    cout << sEmpty.GetDissipator() << endl;
    cout << sUp.GetDissipator() << endl;
    cout << sDown.GetDissipator() << endl;
    return 0;
}
