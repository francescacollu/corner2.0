//
//  testSystem.cpp
//  
//
//  Created by Francesca Collu on 11/11/2018.
//

#include "System.hpp"
#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    TApplication* app = new TApplication("app", 0, 0);
    TCanvas* c = new TCanvas;

    int M = atoi(argv[1]);
    float Jx = atof(argv[2]);
    float Jy = atof(argv[3]);
    float Jz = atof(argv[4]);
    int cod = atoi(argv[5]); //int needed to mark this output

    string s = "1Up_1Down8sM";
    s.append(argv[1]);
    s.append("JxJyJz");
    s.append(argv[5]);
    s.append(".txt");

    const char *output = s.c_str();

    System sy;

    sy.Add(Site(ZUp));
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site(ZDown));
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

    TGraph *p = new TGraph("1Up_1Down8sM60JxJyJz1051.txt");
    p->SetMarkerStyle(20);
    p->SetMarkerColor(kBlack);
    p->SetLineWidth(2);
    p->Draw("ALP");

    c->Update();
    app->Run();

    return 0;
}
