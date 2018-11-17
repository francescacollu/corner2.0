//
//  testSystem.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "System.hpp"
#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TApplication.h>

using namespace std;

int main()
{
    TApplication* app = new TApplication("app", 0, 0);
    TCanvas* c = new TCanvas;

    std::cout << "CIAAAAO" << std::endl;
    TH1F* h = new TH1F;
    h->Fill(0.5);
    h->Draw();

    arma::cx_mat m;

    c->Update();

    System sy;
    
    sy.Add(Site(Site::DissipatorType::ZUp));
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());

    // sy.Add(Site(Site::DissipatorType::ZDown));
    // sy.Add(Site(Site::DissipatorType::ZDown));
    // sy.Add(Site(Site::DissipatorType::ZDown));
    // sy.Add(Site(Site::DissipatorType::ZDown));
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());


    sy.SetCornerSize(40);

    sy.Simulate();

    sy.GetExpValue();

    app->Run();

    return 0;
}
