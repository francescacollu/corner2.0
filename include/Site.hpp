//
//  Site.hpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#ifndef Site_hpp
#define Site_hpp

enum DissipatorType{Empty, ZUp, ZDown};

class Site
{
public:
    DissipatorType GetDissipator() const;
    Site(DissipatorType);
    Site();
    
private:
    DissipatorType Dissipator;
};

#endif /* Site_hpp */
