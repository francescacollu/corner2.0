//
//  Site.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "Site.hpp"

Site::Site(DissipatorType d)
{
    Dissipator = d;
}

Site::Site()
{
    Dissipator = Empty;
}

DissipatorType Site::GetDissipator() const
{
    return Dissipator;
}
