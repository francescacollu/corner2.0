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

    string s = "../out/SpinProfile_s8_M";
    s.append(argv[1]);
    s.append("JxJyJz");
    s.append(argv[5]);
    s.append(".txt");

    const char* output = s.c_str();

    ofstream myfile(output);

    cout << __LINE__ << endl;
    System sy;

    ///////////////////////////////////
    //////////////////////////////////
    //////Convergence study////////////
    // for(int M = 4; M <= 30; M++)
    // {
    //    System sy;

    //    sy.Add(Site(ZUp));
    //    sy.Add(Site());
    //    sy.Add(Site());
    //    sy.Add(Site());
    //    sy.Add(Site());
    //    sy.Add(Site());
    //    sy.Add(Site());
    //    sy.Add(Site(ZDown));

    //    sy.SetCouplingConstants(Jx, Jy, Jz);
    //    cout << "\nM = " << M << "\n";
    //    sy.SetCornerSize(M);
    //    sy.Simulate();
    //    myfile << M << "\t" << sy.Convergence(M) << endl;
    // }
    ///////////////////////////////////
    //////////////////////////////////

    sy.Add(Site(ZUp));
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site(ZDown));

    cout << __LINE__ << endl;
    sy.SetCouplingConstants(Jx, Jy, Jz);
    sy.SetCornerSize(M);
    cout << __LINE__ << endl;
    sy.Simulate();
    cout << __LINE__ << endl;
    sy.GetExpValue(output);
    //sy.Get2PCorrelationFunction(output, 0);


    TGraph *p = new TGraph(output);
    p->SetMarkerStyle(20);
    p->SetMarkerColor(kBlue);
    p->SetLineWidth(2);
    p->Draw("AP");
    p->SetTitle(output);
    c->SetGrid();
    //p->GetXaxis()->SetTitle("M");
    //p->GetYaxis()->SetTitle("<Sz>");

    c->Update();
    c->SaveAs("../../RESULTS/CSR/LM/SpinProfile_s8_M60_1Up_1DownJxJyJz1050.pdf");
    app->Run();

    return 0;
}
