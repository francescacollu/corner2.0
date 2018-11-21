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
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    TApplication* app = new TApplication("app", 0, 0);
    TCanvas* c = new TCanvas;

    TH1F* h = new TH1F;
    h->Fill(0.5);
    h->Draw();

    arma::cx_mat m;

    c->Update();

    int M = atoi(argv[1]);
    float Jx = atof(argv[2]);
    float Jy = atof(argv[3]);
    float Jz = atof(argv[4]);
    int cod = atoi(argv[5]); //int needed to mark this output

    string s = "4sM";
    s.append(argv[1]);
    s.append("JxJyJz");
    s.append(argv[5]);
    s.append(".txt");

    const char *output = s.c_str();

    System sy;

    sy.Add(Site(ZDown));
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site()); 
    // sy.Add(Site());

    sy.SetCouplingConstants(Jx, Jy, Jz);
    sy.SetCornerSize(M);

    sy.Simulate();

    sy.GetExpValue(output);

    app->Run();

    return 0;
}
