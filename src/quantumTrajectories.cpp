#include <armadillo>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    string s = "../../RESULTS/QT/CorrFunc1_s8J105";
    s.append(argv[1]);
    s.append(".txt");
    const char* filename = s.c_str();
    ofstream myfileCorrFunc1(filename);

    string s2 = "../../RESULTS/QT/LM_s8J105";
    s2.append(argv[1]);
    s2.append(".txt");
    const char* filename2 = s2.c_str();
    ofstream myfileMagnetiz(filename2);

    string s3 = "../../RESULTS/QT/spinCurr_s8J105";
    s3.append(argv[1]);
    s3.append(".txt");
    const char* filename3 = s3.c_str();
    ofstream myfileSpinCurr(filename3);

    TApplication* app = new TApplication("app", 0, 0);
    TH1F *realZ0 = new TH1F("realZ[0]", "realZ[0]", 100, -7, 7);
    TH1F *realZ1 = new TH1F("realZ[1]", "realZ[1]", 100, -7, 7);
    TH1F *realZ2 = new TH1F("realZ[2]", "realZ[2]", 100, -7, 7);
    TH1F *realZ3 = new TH1F("realZ[3]", "realZ[3]", 100, -7, 7);
    TH1F *realZ4 = new TH1F("realZ[4]", "realZ[4]", 100, -7, 7);
    TH1F *realZ5 = new TH1F("realZ[5]", "realZ[5]", 100, -7, 7);
    TH1F *realZ6 = new TH1F("realZ[6]", "realZ[6]", 100, -7, 7);
    TH1F *realZ7 = new TH1F("realZ[7]", "realZ[7]", 100, -7, 7);

    TH1F *corrFunc00 = new TH1F("corrFunc0[0]", "corrFunc0[0]", 100, -7, 7);
    TH1F *corrFunc01 = new TH1F("corrFunc0[1]", "corrFunc0[1]", 100, -7, 7);
    TH1F *corrFunc02 = new TH1F("corrFunc0[2]", "corrFunc0[2]", 100, -7, 7);
    TH1F *corrFunc03 = new TH1F("corrFunc0[3]", "corrFunc0[3]", 100, -7, 7);
    TH1F *corrFunc04 = new TH1F("corrFunc0[4]", "corrFunc0[4]", 100, -7, 7);
    TH1F *corrFunc05 = new TH1F("corrFunc0[5]", "corrFunc0[5]", 100, -7, 7);
    TH1F *corrFunc06 = new TH1F("corrFunc0[6]", "corrFunc0[6]", 100, -7, 7);
    TH1F *corrFunc07 = new TH1F("corrFunc0[7]", "corrFunc0[7]", 100, -7, 7);

    TH1F *spinCurr01 = new TH1F("spinCurr[0]", "spinCurr[0]", 100, -1, 1);
    TH1F *spinCurr12 = new TH1F("spinCurr[1]", "spinCurr[1]", 100, -1, 1);
    TH1F *spinCurr23 = new TH1F("spinCurr[2]", "spinCurr[2]", 100, -1, 1);
    TH1F *spinCurr34 = new TH1F("spinCurr[3]", "spinCurr[3]", 100, -1, 1);
    TH1F *spinCurr45 = new TH1F("spinCurr[4]", "spinCurr[4]", 100, -1, 1);
    TH1F *spinCurr56 = new TH1F("spinCurr[5]", "spinCurr[5]", 100, -1, 1);
    TH1F *spinCurr67 = new TH1F("spinCurr[6]", "spinCurr[6]", 100, -1, 1);
    
    double Jx = 1., Jy = 0.5;
    float Jz = atof(argv[1]);
    double T = 100.;
    double dt = 0.1;

    double n = 0.;
    int Nmax = round(T/dt);

    cx_mat Sx, Sy, Sz, I;
    std::vector<cx_mat> sigmaZ(8);
    std::vector<cx_mat> sigmaX(8);
    std::vector<cx_mat> sigmaY(8);

    Sx << cx_double(0., 0.) << cx_double(1., 0.) << endr
        << cx_double(1., 0.) << cx_double(0., 0.) << endr;
    
    Sy << cx_double(0., 0.) << cx_double(0., -1.) << endr
        << cx_double(0., 1.) << cx_double(0., 0.) << endr;
    
    Sz << cx_double(1., 0.) << cx_double(0., 0.) << endr
        << cx_double(0, 0.) << cx_double(-1., 0.) << endr;

    I << cx_double(1., 0.) << cx_double(0., 0.) << endr
        << cx_double(0., 0.) << cx_double(1., 0.) << endr;

    cx_mat II;
    II = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));

    sigmaZ[0] = kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaZ[1] = kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaZ[2] = kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaZ[3] = kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, I)))))));
    sigmaZ[4] = kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, I)))))));
    sigmaZ[5] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, I)))))));
    sigmaZ[6] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, I)))))));
    sigmaZ[7] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, Sz)))))));

    sigmaX[0] = kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaX[1] = kron(I, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaX[2] = kron(I, kron(I, kron(Sx, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaX[3] = kron(I, kron(I, kron(I, kron(Sx, kron(I, kron(I, kron(I, I)))))));
    sigmaX[4] = kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(I, kron(I, I)))))));
    sigmaX[5] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(I, I)))))));
    sigmaX[6] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, I)))))));
    sigmaX[7] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, Sx)))))));

    sigmaY[0] = kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaY[1] = kron(I, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaY[2] = kron(I, kron(I, kron(Sy, kron(I, kron(I, kron(I, kron(I, I)))))));
    sigmaY[3] = kron(I, kron(I, kron(I, kron(Sy, kron(I, kron(I, kron(I, I)))))));
    sigmaY[4] = kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(I, kron(I, I)))))));
    sigmaY[5] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(I, I)))))));
    sigmaY[6] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, I)))))));
    sigmaY[7] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, Sy)))))));

    cx_mat Cdown, Cup, K;

    Cdown << cx_double(0., 0.) << cx_double(0., 0.) << endr
        << cx_double(1., 0.) << cx_double(0., 0.) << endr;

    Cup << cx_double(0., 0.) << cx_double(1., 0.) << endr
        << cx_double(0., 0.) << cx_double(0., 0.) << endr;

    Cup = kron(Cup, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))));
    Cdown = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, Cdown)))))));

    K = -0.5 * (arma::trans(Cup)*Cup + arma::trans(Cdown)*Cdown);

    cx_mat H;


    H = Jx*( kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, Sx))))))) ) + 
          Jy*( kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, Sy))))))) ) +
          Jz*( kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, I))))))) + 
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, Sz))))))) );

    cx_mat Heff = H + cx_double(0., 1.)*K;
    cx_mat E = expmat(cx_double(0.,-1.)*Heff*dt);

    cx_mat eigvecH;
    cx_vec eigvalH, sorted_eigvalH;

    eig_gen(eigvalH, eigvecH, H);
    uvec index_eigval = stable_sort_index(eigvalH);
    eigvecH = eigvecH.cols(index_eigval);
    sorted_eigvalH = sort(eigvalH);

    cx_vec phi0, phi;
    cx_rowvec phi0T;
    std::vector<cx_vec> PHI;

    std::vector<cx_double> SigmaZ0;
    std::vector<cx_double> SigmaZ1;
    std::vector<cx_double> SigmaZ2;
    std::vector<cx_double> SigmaZ3;
    std::vector<cx_double> SigmaZ4;
    std::vector<cx_double> SigmaZ5;
    std::vector<cx_double> SigmaZ6;
    std::vector<cx_double> SigmaZ7;

    std::vector<cx_double> corrFunc0;
    std::vector<cx_double> corrFunc1;
    std::vector<cx_double> corrFunc2;
    std::vector<cx_double> corrFunc3;
    std::vector<cx_double> corrFunc4;
    std::vector<cx_double> corrFunc5;
    std::vector<cx_double> corrFunc6;
    std::vector<cx_double> corrFunc7;

    std::vector<cx_double> spinCurr0;
    std::vector<cx_double> spinCurr1;
    std::vector<cx_double> spinCurr2;
    std::vector<cx_double> spinCurr3;
    std::vector<cx_double> spinCurr4;
    std::vector<cx_double> spinCurr5;
    std::vector<cx_double> spinCurr6;
    
    phi0 = eigvecH.col(0);

    cx_mat DPUp, DPDown;
    double dp, dpUp, dpDown;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0., 1.);
    
    int N=0;
    while(N<1000)
    {   
        N++;
        cout << "N = " << N << endl;
        for(int n=0; n<Nmax; n++)
        {
            // cout << endl << n << endl;
            DPUp = (arma::trans(Cup*phi0))*(Cup*phi0)*dt;
            DPDown = (arma::trans(Cdown*phi0))*(Cdown*phi0)*dt;

            dpUp = real(DPUp.eval().at(0,0));
            dpDown = real(DPDown.eval().at(0,0));
            dp = dpUp + dpDown;

            // cout << "dpUp = " << dpUp << endl;
            // cout << "dpDown = " << dpDown << endl;
            double eps = distribution(generator);

            if(eps>0. && eps<=dpUp)
            {
                phi = Cup*phi0/sqrt(cdot(Cup*phi0, Cup*phi0));
                //std::cout << "eps = " << eps << endl << "dpUp = " << dpUp << endl;
            }
            else if(eps>dpUp && eps<=dp)
            {
                phi = Cdown*phi0/sqrt(cdot(Cdown*phi0, Cdown*phi0));
                //std::cout << "eps = " << eps << endl << "dpDown = " << dpDown << endl;
            }
            else if(eps>dp && eps<=1.)
            {
                phi = E*phi0/sqrt(1.-dp);
                phi = phi/sqrt((trans(phi)*phi).eval().at(0,0));
                //std::cout << "eps = " << eps << endl << "dp = " << dp << endl;
            }
            phi = phi/sqrt((trans(phi)*phi).eval().at(0,0));
            phi0 = phi;
        }

        PHI.push_back(phi);

        SigmaZ0.push_back((arma::trans(PHI.back())*sigmaZ[0]*PHI.back()).eval().at(0,0));
        realZ0->Fill(real(SigmaZ0.back()));

        SigmaZ1.push_back((arma::trans(PHI.back())*sigmaZ[1]*PHI.back()).eval().at(0,0));
        realZ1->Fill(real(SigmaZ1.back()));

        SigmaZ2.push_back((arma::trans(PHI.back())*sigmaZ[2]*PHI.back()).eval().at(0,0));
        realZ2->Fill(real(SigmaZ2.back()));

        SigmaZ3.push_back((arma::trans(PHI.back())*sigmaZ[3]*PHI.back()).eval().at(0,0));
        realZ3->Fill(real(SigmaZ3.back()));

        SigmaZ4.push_back((arma::trans(PHI.back())*sigmaZ[4]*PHI.back()).eval().at(0,0));
        realZ4->Fill(real(SigmaZ4.back()));

        SigmaZ5.push_back((arma::trans(PHI.back())*sigmaZ[5]*PHI.back()).eval().at(0,0));
        realZ5->Fill(real(SigmaZ5.back()));

        SigmaZ6.push_back((arma::trans(PHI.back())*sigmaZ[6]*PHI.back()).eval().at(0,0));
        realZ6->Fill(real(SigmaZ6.back()));

        SigmaZ7.push_back((arma::trans(PHI.back())*sigmaZ[7]*PHI.back()).eval().at(0,0));
        realZ7->Fill(real(SigmaZ7.back()));

        corrFunc0.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[0]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[0]*PHI.back())).eval().at(0,0));
        corrFunc00->Fill(real(corrFunc0.back()));

        corrFunc1.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[1]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[1]*PHI.back())).eval().at(0,0));
        corrFunc01->Fill(real(corrFunc1.back()));

        corrFunc2.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[2]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[2]*PHI.back())).eval().at(0,0));
        corrFunc02->Fill(real(corrFunc2.back()));

        corrFunc3.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[3]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[3]*PHI.back())).eval().at(0,0));
        corrFunc03->Fill(real(corrFunc3.back()));

        corrFunc4.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[4]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[4]*PHI.back())).eval().at(0,0));
        corrFunc04->Fill(real(corrFunc4.back()));

        corrFunc5.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[5]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[5]*PHI.back())).eval().at(0,0));
        corrFunc05->Fill(real(corrFunc5.back()));

        corrFunc6.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[6]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[6]*PHI.back())).eval().at(0,0));
        corrFunc06->Fill(real(corrFunc6.back()));

        corrFunc7.push_back((trans(PHI.back())*sigmaZ[0]*sigmaZ[7]*PHI.back()).eval().at(0,0));
        //-((trans(PHI.back())*sigmaZ[0]*PHI.back())*(trans(PHI.back())*sigmaZ[7]*PHI.back())).eval().at(0,0));
        corrFunc07->Fill(real(corrFunc7.back()));


        spinCurr0.push_back((trans(PHI.back())*(sigmaX[0]*sigmaY[1]-sigmaY[0]*sigmaX[1])*PHI.back()).eval().at(0,0));
        spinCurr01->Fill(real(spinCurr0.back()));

        spinCurr1.push_back((trans(PHI.back())*(sigmaX[1]*sigmaY[2]-sigmaY[1]*sigmaX[2])*PHI.back()).eval().at(0,0));
        spinCurr12->Fill(real(spinCurr1.back()));

        spinCurr2.push_back((trans(PHI.back())*(sigmaX[2]*sigmaY[3]-sigmaY[2]*sigmaX[3])*PHI.back()).eval().at(0,0));
        spinCurr23->Fill(real(spinCurr2.back()));

        spinCurr3.push_back((trans(PHI.back())*(sigmaX[3]*sigmaY[4]-sigmaY[3]*sigmaX[4])*PHI.back()).eval().at(0,0));
        spinCurr34->Fill(real(spinCurr3.back()));

        spinCurr4.push_back((trans(PHI.back())*(sigmaX[4]*sigmaY[5]-sigmaY[4]*sigmaX[5])*PHI.back()).eval().at(0,0));
        spinCurr45->Fill(real(spinCurr4.back()));

        spinCurr5.push_back((trans(PHI.back())*(sigmaX[5]*sigmaY[6]-sigmaY[5]*sigmaX[6])*PHI.back()).eval().at(0,0));
        spinCurr56->Fill(real(spinCurr5.back()));

        spinCurr6.push_back((trans(PHI.back())*(sigmaX[6]*sigmaY[7]-sigmaY[6]*sigmaX[7])*PHI.back()).eval().at(0,0));
        spinCurr67->Fill(real(spinCurr6.back()));
    }
    std::cout << "PHI.size = " << PHI.size() << endl;
    std::cout << "N = " << N << endl;

    myfileMagnetiz << 1 << '\t' << realZ0->GetMean() << endl;
    myfileMagnetiz << 2 << '\t' << realZ1->GetMean() << endl;
    myfileMagnetiz << 3 << '\t' << realZ2->GetMean() << endl;
    myfileMagnetiz << 4 << '\t' << realZ3->GetMean() << endl;
    myfileMagnetiz << 5 << '\t' << realZ4->GetMean() << endl;
    myfileMagnetiz << 6 << '\t' << realZ5->GetMean() << endl;
    myfileMagnetiz << 7 << '\t' << realZ6->GetMean() << endl;
    myfileMagnetiz << 8 << '\t' << realZ7->GetMean() << endl;

    myfileCorrFunc1 << 1 << '\t' << corrFunc00->GetMean() << endl;
    myfileCorrFunc1 << 2 << '\t' << corrFunc01->GetMean() << endl;
    myfileCorrFunc1 << 3 << '\t' << corrFunc02->GetMean() << endl;
    myfileCorrFunc1 << 4 << '\t' << corrFunc03->GetMean() << endl;
    myfileCorrFunc1 << 5 << '\t' << corrFunc04->GetMean() << endl;
    myfileCorrFunc1 << 6 << '\t' << corrFunc05->GetMean() << endl;
    myfileCorrFunc1 << 7 << '\t' << corrFunc06->GetMean() << endl;
    myfileCorrFunc1 << 8 << '\t' << corrFunc07->GetMean() << endl;

    myfileSpinCurr << 1 << '\t' << spinCurr01->GetMean() << endl;
    myfileSpinCurr << 2 << '\t' << spinCurr12->GetMean() << endl;
    myfileSpinCurr << 3 << '\t' << spinCurr23->GetMean() << endl;
    myfileSpinCurr << 4 << '\t' << spinCurr34->GetMean() << endl;
    myfileSpinCurr << 5 << '\t' << spinCurr45->GetMean() << endl;
    myfileSpinCurr << 6 << '\t' << spinCurr56->GetMean() << endl;
    myfileSpinCurr << 7 << '\t' << spinCurr67->GetMean() << endl;

    TCanvas* can1 = new TCanvas;
    spinCurr01->Draw();
    can1->Update();

    TCanvas* can2 = new TCanvas;
    spinCurr12->Draw();
    can2->Update();

    TCanvas* can3 = new TCanvas;
    spinCurr23->Draw();
    can3->Update();

    TCanvas* can4 = new TCanvas;
    spinCurr34->Draw();
    can4->Update();

    TCanvas* can5 = new TCanvas;
    spinCurr45->Draw();
    can5->Update();

    TCanvas* can6 = new TCanvas;
    spinCurr56->Draw();
    can6->Update();

    TCanvas* can7 = new TCanvas;
    spinCurr67->Draw();
    can7->Update();

    TCanvas* c = new TCanvas;
    TGraph* p = new TGraph(filename);
    p->SetMarkerStyle(20);
    p-> Draw("ALP");
    c -> Update();

    TCanvas* c2 = new TCanvas;
    TGraph* p2 = new TGraph(filename2);
    p2->SetMarkerStyle(20);
    p2-> Draw("ALP");
    c2 -> Update();

    TCanvas* c3 = new TCanvas;
    TGraph* p3 = new TGraph(filename3);
    p3->SetMarkerStyle(20);
    p3-> Draw("ALP");
    c3 -> Update();

    app->Run();

    delete realZ0;
    delete realZ1;
    delete realZ2;
    delete realZ3;
    delete realZ4;
    delete realZ5;
    delete realZ6;
    delete realZ7;

    delete corrFunc00;
    delete corrFunc01;
    delete corrFunc02;
    delete corrFunc03;
    delete corrFunc04;
    delete corrFunc05;
    delete corrFunc06;
    delete corrFunc07;

    delete spinCurr01;
    delete spinCurr12;
    delete spinCurr23;
    delete spinCurr34;
    delete spinCurr45;
    delete spinCurr56;
    delete spinCurr67;

    delete can1;
    delete can2;
    delete can3;
    delete can4;
    delete can5;
    delete can6;
    delete can7;

    delete c;
    delete c2;
    delete c3;
    delete app;

    return 0;
}