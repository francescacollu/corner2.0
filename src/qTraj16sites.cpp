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

int main()
{
    TApplication* app = new TApplication("app", 0, 0);
    TH1F *realZ0 = new TH1F("realZ[0]", "realZ[0]", 100, -1., 1.);
    TH1F *realZ1 = new TH1F("realZ[1]", "realZ[1]", 100, -1., 1.);
    TH1F *realZ2 = new TH1F("realZ[2]", "realZ[2]", 100, -1., 1.);
    TH1F *realZ3 = new TH1F("realZ[3]", "realZ[3]", 100, -1., 1.);
    TH1F *realZ4 = new TH1F("realZ[4]", "realZ[4]", 100, -1., 1.);
    TH1F *realZ5 = new TH1F("realZ[5]", "realZ[5]", 100, -1., 1.);
    TH1F *realZ6 = new TH1F("realZ[6]", "realZ[6]", 100, -1., 1.);
    TH1F *realZ7 = new TH1F("realZ[7]", "realZ[7]", 100, -1., 1.);
    TH1F *realZ8 = new TH1F("realZ[8]", "realZ[8]", 100, -1., 1.);
    TH1F *realZ9 = new TH1F("realZ[9]", "realZ[9]", 100, -1., 1.);
    TH1F *realZ10 = new TH1F("realZ[10]", "realZ[10]", 100, -1., 1.);
    TH1F *realZ11 = new TH1F("realZ[11]", "realZ[11]", 100, -1., 1.);
    TH1F *realZ12 = new TH1F("realZ[12]", "realZ[12]", 100, -1., 1.);
    TH1F *realZ13 = new TH1F("realZ[13]", "realZ[13]", 100, -1., 1.);
    TH1F *realZ14 = new TH1F("realZ[14]", "realZ[14]", 100, -1., 1.);
    TH1F *realZ15 = new TH1F("realZ[15]", "realZ[15]", 100, -1., 1.);
    
    double Jx = 1., Jy = 0.5, Jz = 1.;
    double T = 100.;
    double dt = 0.1;

    double n = 0.;
    int Nmax = round(T/dt);

    cx_mat Sx, Sy, Sz, I;
    std::vector<cx_mat> sigmaZ(16);

    Sx << cx_double(0., 0.) << cx_double(0.5, 0.) << endr
        << cx_double(0.5, 0.) << cx_double(0., 0.) << endr;
    
    Sy << cx_double(0., 0.) << cx_double(0., -0.5) << endr
        << cx_double(0., 0.5) << cx_double(0., 0.) << endr;
    
    Sz << cx_double(0.5, 0.) << cx_double(0., 0.) << endr
        << cx_double(0, 0.) << cx_double(-0.5, 0.) << endr;

    I << cx_double(1., 0.) << cx_double(0., 0.) << endr
        << cx_double(0., 0.) << cx_double(1., 0.) << endr;
    
    // sigmaZ[0] = kron(Sz, kron(I, kron(I, I)));
    // sigmaZ[1] = kron(I, kron(Sz, kron(I, I)));
    // sigmaZ[2] = kron(I, kron(I, kron(Sz, I)));
    // sigmaZ[3] = kron(I, kron(I, kron(I, Sz)));

    cout << "16x16-identity-matrix...\n";
    cx_mat II;
    II = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));

    cout << "16x16-sigmaZ-matrices...\n";
    sigmaZ[0] = kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[1] = kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[2] = kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[3] = kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[4] = kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[5] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[6] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[7] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[8] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[9] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[10] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[11] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, kron(I, I)))))))))))))));
    sigmaZ[12] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, kron(I, I)))))))))))))));
    sigmaZ[13] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(I, I)))))))))))))));
    sigmaZ[14] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, I)))))))))))))));
    sigmaZ[15] = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, Sz)))))))))))))));

    cx_mat Cdown, Cup, K;

    Cdown << cx_double(0., 0.) << cx_double(0., 0.) << endr
        << cx_double(1., 0.) << cx_double(0., 0.) << endr;

    Cup << cx_double(0., 0.) << cx_double(1., 0.) << endr
        << cx_double(0., 0.) << cx_double(0., 0.) << endr;

    // Cup = kron(Cup, kron(I, kron(I, I)));
    // Cdown = kron(I, kron(I, kron(I, Cdown)));
    Cup = kron(Cup, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I)))))))))))))));
    Cdown = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, Cdown)))))))))))))));

    K = -0.5 * (arma::trans(Cup)*Cup + arma::trans(Cdown)*Cdown);

    cx_mat H;
    // H = - Jx* (kron(Sx, kron(Sx, kron(I, I))) + kron(I, kron(Sx, kron(Sx, I))) +//
    //      kron(I, kron(I, kron(Sx, Sx)))) - Jy * (kron(Sy, kron(Sy, kron(I, I))) +//
    //      kron(I, kron(Sy, kron(Sy, I))) + kron(I, kron(I, kron(Sy, Sy)))) -//
    //      Jz * (kron(Sz, kron(Sz, kron(I, I))) +  kron(I, kron(Sz, kron(Sz, I))) +//
    //      kron(I, kron(I, kron(Sz, Sz))));


    H = - Jx*( kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, kron(Sx, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sx, Sx))))))))))))))) ) -
          Jy*( kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, kron(Sy, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sy, Sy))))))))))))))) ) - 
          Jz*( kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, kron(I, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, kron(Sz, I))))))))))))))) +
               kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(Sz, Sz))))))))))))))) );

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
    std::vector<cx_double> SigmaZ8;
    std::vector<cx_double> SigmaZ9;
    std::vector<cx_double> SigmaZ10;
    std::vector<cx_double> SigmaZ11;
    std::vector<cx_double> SigmaZ12;
    std::vector<cx_double> SigmaZ13;
    std::vector<cx_double> SigmaZ14;
    std::vector<cx_double> SigmaZ15;
    
    phi0 = eigvecH.col(0);

    cx_mat DPUp, DPDown;
    double dp, dpUp, dpDown;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0., 1.);
    
    int N=0;
    while(N<10)
    {   
        N++;
        cout << "N = " << N << endl;
        for(int n=0; n<Nmax; n++)
        {
            cout << endl << n << endl;
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

        SigmaZ8.push_back((arma::trans(PHI.back())*sigmaZ[8]*PHI.back()).eval().at(0,0));
        realZ8->Fill(real(SigmaZ8.back()));

        SigmaZ9.push_back((arma::trans(PHI.back())*sigmaZ[9]*PHI.back()).eval().at(0,0));
        realZ9->Fill(real(SigmaZ9.back()));

        SigmaZ10.push_back((arma::trans(PHI.back())*sigmaZ[10]*PHI.back()).eval().at(0,0));
        realZ10->Fill(real(SigmaZ10.back()));

        SigmaZ11.push_back((arma::trans(PHI.back())*sigmaZ[11]*PHI.back()).eval().at(0,0));
        realZ11->Fill(real(SigmaZ11.back()));

        SigmaZ12.push_back((arma::trans(PHI.back())*sigmaZ[12]*PHI.back()).eval().at(0,0));
        realZ12->Fill(real(SigmaZ12.back()));

        SigmaZ13.push_back((arma::trans(PHI.back())*sigmaZ[13]*PHI.back()).eval().at(0,0));
        realZ13->Fill(real(SigmaZ13.back()));

        SigmaZ14.push_back((arma::trans(PHI.back())*sigmaZ[14]*PHI.back()).eval().at(0,0));
        realZ14->Fill(real(SigmaZ14.back()));

        SigmaZ15.push_back((arma::trans(PHI.back())*sigmaZ[15]*PHI.back()).eval().at(0,0));
        realZ15->Fill(real(SigmaZ15.back()));

    }
    std::cout << "PHI.size = " << PHI.size() << endl;
    std::cout << "N = " << N << endl;

    TCanvas* c0 = new TCanvas;
    realZ0->Draw();
    c0->Update();

    TCanvas* c1 = new TCanvas;
    realZ1->Draw();
    c1->Update();

    TCanvas* c2 = new TCanvas;
    realZ2->Draw();
    c2->Update();

    TCanvas* c3 = new TCanvas;
    realZ3->Draw();
    c3->Update();

    TCanvas* c4 = new TCanvas;
    realZ4->Draw();
    c4->Update();

    TCanvas* c5 = new TCanvas;
    realZ5->Draw();
    c5->Update();

    TCanvas* c6 = new TCanvas;
    realZ6->Draw();
    c6->Update();

    TCanvas* c7 = new TCanvas;
    realZ7->Draw();
    c7->Update();

    TCanvas* c8 = new TCanvas;
    realZ8->Draw();
    c8->Update();

    TCanvas* c9 = new TCanvas;
    realZ9->Draw();
    c9->Update();

    TCanvas* c10 = new TCanvas;
    realZ10->Draw();
    c10->Update();

    TCanvas* c11 = new TCanvas;
    realZ11->Draw();
    c11->Update();

    TCanvas* c12 = new TCanvas;
    realZ12->Draw();
    c12->Update();

    TCanvas* c13 = new TCanvas;
    realZ13->Draw();
    c13->Update();

    TCanvas* c14 = new TCanvas;
    realZ14->Draw();
    c14->Update();

    TCanvas* c15 = new TCanvas;
    realZ15->Draw();
    c15->Update();


    c0->SaveAs("QTsigmaZ[0].pdf");
    c1->SaveAs("QTsigmaZ[1].pdf");
    c2->SaveAs("QTsigmaZ[2].pdf");
    c3->SaveAs("QTsigmaZ[3].pdf");
    c4->SaveAs("QTsigmaZ[4].pdf");
    c5->SaveAs("QTsigmaZ[5].pdf");
    c6->SaveAs("QTsigmaZ[6].pdf");
    c7->SaveAs("QTsigmaZ[7].pdf");
    c8->SaveAs("QTsigmaZ[8].pdf");
    c9->SaveAs("QTsigmaZ[9].pdf");
    c10->SaveAs("QTsigmaZ[10].pdf");
    c11->SaveAs("QTsigmaZ[11].pdf");
    c12->SaveAs("QTsigmaZ[12].pdf");
    c13->SaveAs("QTsigmaZ[13].pdf");
    c14->SaveAs("QTsigmaZ[14].pdf");
    c15->SaveAs("QTsigmaZ[15].pdf");
    app->Run();

    delete realZ0;
    delete realZ1;
    delete realZ2;
    delete realZ3;
    delete realZ4;
    delete realZ5;
    delete realZ6;
    delete realZ7;
    delete realZ8;
    delete realZ9;
    delete realZ10;
    delete realZ11;
    delete realZ12;
    delete realZ13;
    delete realZ14;
    delete realZ15;

    delete c0;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete c6;
    delete c7;
    delete c8;
    delete c9;
    delete c10;
    delete c11;
    delete c12;
    delete c13;
    delete c14;
    delete c15;
    delete app;


    return 0;
}