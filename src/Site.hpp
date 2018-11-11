//
//  Site.hpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#ifndef Site_hpp
#define Site_hpp

class Site
{
public:
    enum DissipatorType{Empty, ZUp, ZDown};
    DissipatorType GetDissipator() const;
    Site(DissipatorType);
    Site();
    
private:
    DissipatorType Dissipator;
};

#endif /* Site_hpp */
