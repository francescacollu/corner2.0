#include "System.hpp"
#include <iostream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    TApplication* app = new TApplication("app", 0, 0);

    int M = atoi(argv[1]);
    float Jx = atof(argv[2]);
    float Jy = atof(argv[3]);
    float Jz = atof(argv[4]);

    string s = "../../RESULTS/CSR/s16_M";
    s.append(argv[1]);
    s.append("J");
    s.append(argv[2]);
    s.append(argv[3]);
    s.append(argv[4]);

    // const char* out = s.c_str();

    // string stxt = s;
    // stxt.append(".txt");

    // const char* outTxt = stxt.c_str();

    // ofstream myfile(outTxt);

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
    ///////////////////////////////////

    // sy.Add(Site(ZUp));
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    // sy.Add(Site());
    sy.Add(Site(ZUp));
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site());
    sy.Add(Site(ZDown));

    sy.SetCouplingConstants(Jx, Jy, Jz);
    sy.SetCornerSize(M);
    sy.Simulate();

    string sLM = s;
    sLM.append("_LM.txt");
    const char* outLM = sLM.c_str();
    sy.GetExpValue(outLM);

    string sCorrFunc = s;
    sCorrFunc.append("_CorrFunc.txt");
    const char* outCorrFunc = sCorrFunc.c_str();
    sy.Get2PCorrelationFunction(outCorrFunc);

    string sSpinCurr = s;
    sSpinCurr.append("_SpinCurr.txt"); 
    const char* outSpinCurr = sSpinCurr.c_str();
    sy.GetSpinCurrent(outSpinCurr);

    TCanvas* cLM = new TCanvas;
    TGraph *pLM = new TGraph(outLM);
    pLM->SetMarkerStyle(20);
    pLM->SetMarkerColor(kBlue);
    pLM->SetLineWidth(1);
    pLM->Draw("ALP");
    pLM->SetTitle(outLM);
    pLM->GetXaxis()->SetTitle("site_index");
    pLM->GetYaxis()->SetTitle("<Sz>");
    cLM->SetGrid();
    cLM->Update();

    TCanvas* cSC = new TCanvas;
    TGraph *pSC = new TGraph(outSpinCurr);
    pSC->SetMarkerStyle(20);
    pSC->SetMarkerColor(kBlue);
    pSC->SetLineWidth(1);
    pSC->Draw("ALP");
    pSC->SetTitle(outSpinCurr);
    pSC->GetXaxis()->SetTitle("site_index");
    pSC->GetYaxis()->SetTitle("<J>");
    cSC->SetGrid();
    cSC->Update();

    string sLMpdf = s;
    sLMpdf.append("_LM.pdf");
    const char* outLMPdf = sLMpdf.c_str();
    cLM->SaveAs(outLMPdf);

    string sSCpdf = s;
    sSCpdf.append("_SpinCurr.pdf");
    const char* outSCPdf = sSCpdf.c_str();
    cLM->SaveAs(outSCPdf);

    app->Run();

    delete cLM;

    return 0;
}
