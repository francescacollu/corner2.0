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

    //int M = atoi(argv[1]);
    float Jx = atof(argv[1]);
    float Jy = atof(argv[2]);
    float Jz = atof(argv[3]);
    int cod = atoi(argv[4]); //int needed to mark this output

    string s = "../out/prova";
    s.append(argv[1]);
    s.append("JxJyJz");
    s.append(argv[4]);
    s.append(".txt");

    const char *output = s.c_str();

    ofstream myfile(output);

    System sy;

    ///////////////////////////////////
    //////////////////////////////////
    //////Convergence study////////////
    for(int M = 4; M <= 30; M++)
    {
       System sy;

       sy.Add(Site(ZUp));
       sy.Add(Site());
       sy.Add(Site());
       sy.Add(Site());
       sy.Add(Site());
       sy.Add(Site());
       sy.Add(Site());
       sy.Add(Site(ZDown));

       sy.SetCouplingConstants(Jx, Jy, Jz);
       cout << "\nM = " << M << "\n";
       sy.SetCornerSize(M);
       sy.Simulate();
       myfile << M << "\t" << sy.Convergence(M) << endl;
    }
    ///////////////////////////////////
    //////////////////////////////////

    // sy.Add(Site(ZUp));
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site(ZDown));

    // sy.SetCouplingConstants(Jx, Jy, Jz);
    // sy.SetCornerSize(M);
    // sy.Simulate();
    // sy.GetExpValue(output);
    //sy.Get2PCorrelationFunction(output, 0);


    TGraph *p = new TGraph(output);
    //TGraph *p = new TGraph("corrFunc_8sites.txt");
    p->SetMarkerStyle(20);
    p->SetMarkerColor(kBlue);
    p->SetLineWidth(2);
    p->Draw("AP");
    p->SetTitle(output);
    //p->GetXaxis()->SetTitle("M");
    //p->GetYaxis()->SetTitle("<Sz>");

    c->Update();
    c->SaveAs("../../RESULTS/CSR/convergence/CorrFunc_M60_1Up_1Down8sJxJyJz1051.pdf");
    app->Run();

    return 0;
}
