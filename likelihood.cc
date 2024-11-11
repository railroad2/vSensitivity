#include <iostream>
#include "TROOT.h"
#include "TF1.h"

#include "../vOscillation/vClass.hh"

#define DEBUG 0

//TF1* getPDF0()
vInterpolator* getPDF0()
{

    vDetectorCylinder *det = new vDetectorCylinder();
    det->SetRadius(7.5);
    det->SetHeight(15);

    double Lmin = 1.8;
    double Lmax = 20.;

    int npts = 1000;

    vector<double> x, y; 

    for (int i=0; i<npts; i++) {
        double L = (Lmax - Lmin)/npts * i + Lmin;
        x.push_back(L);
        y.push_back(det->GetLDistribution(L, 10, 0, 1001, true)/ 0.747);
        if (DEBUG > 2) cout << L << ' ' << y[i] << endl;
    }

    vInterpolator *interp = new vInterpolator(x, y);

    return interp;
}

vInterpolator* getPDF1(double sin2theta14, double Dm2_41)
{

    vDetectorCylinder *det = new vDetectorCylinder();
    det->SetRadius(7.5);
    det->SetHeight(15);

    double Lmin = 1.8;
    double Lmax = 20.;

    int npts = 1000;

    vSterile* vst = new vSterile();
    vst->LoadStdData();
    vst->Load4StdData();

    double theta14 = TMath::ASin(TMath::Sqrt(sin2theta14));
    double theta24 = TMath::ASin(TMath::Sqrt(0.15));
    double theta34 = 0.;

    vst->Set4radian(theta14, theta24, theta34); // (theta_14, theta_24, theta_34)
    vst->Set4Dm2(Dm2_41);

    vector<double> x, y; 

    for (int i=0; i<npts; i++) {
        double L = (Lmax - Lmin)/npts * i + Lmin;
        double tmp = det->GetLDistribution(L, 10, 0, 1001, true);
        double osc = vst->GetProbability(L, 0.747, "e", "e");
        x.push_back(L);
        y.push_back(tmp * osc);
        if (DEBUG > 2) cout << L << "    " << tmp << "    " << osc << "    " << y[i] << endl;
    }

    vInterpolator *interp = new vInterpolator(x, y);

    return interp;
}

double get_oscillation(double theta_14, double Dm2_41)
{
    return 0;
}

double* readfile(TString fn)
{
    ifstream infile;
    infile.open(fn);
    double L, LoE;

    double arr[10000];

    int nlines = 0;

    while (1) {
        infile >> L >> LoE;
        //cout << LoE << endl;
        if (!infile.good()) 
            break;

        arr[nlines] = LoE;
        nlines++;
    }

    cerr << nlines << " are read." << endl;

    infile.close();

    return &(arr[0]);
}

double logLikelihood_unbinned(double* data, TF1* pdf, double sin2theta14=0, double Dm2_41=0)
{
    double llh = 0;
    double lh = 0;
    int nbin = 100;

    vSterile* vst = new vSterile();
    vst->Load4StdData();

    double theta14 = TMath::ASin(TMath::Sqrt(sin2theta14));
    double theta24 = TMath::ASin(TMath::Sqrt(0.15));
    double theta34 = 0;

    vst->Set4radian(theta14, theta24, theta34); 
    vst->Set4Dm2(Dm2_41);

    for (int i=0; i<nbin; i++) {
        double tmp = pdf->Eval(data[i]/1000.);
        if (tmp <= 0) {
            tmp = 1e-31;
        }
        double osc = vst->GetProbability(data[i]/1000., 1, "e", "e");

        if (DEBUG>2) cout << data[i]/1000. << '\t'  <<  tmp << '\t' << osc << endl;

        llh += TMath::Log(tmp * osc);
    }
    
    return llh;
}

void operation_check()
{
    // getting pdf w/o oscillation
    vInterpolator* interp0 = getPDF0();
    TF1* pdf0 = interp0->GetTF1();

    // getting pdf with oscillation
    double Dm2_41 = 0.01;
    double theta_14 = 0.01;
    vInterpolator* interp1 = getPDF1(theta_14, Dm2_41);
    TF1* pdf1 = interp1->GetTF1();

    pdf0->Draw();
    pdf1->Draw("SAME");
}

void compute_likelihood_org()
{
    //double data[100] = {0,};

    vInterpolator* interp = getPDF0();
    TF1* pdf = interp->GetTF1();

    TString fn = "Cr51_L,LoE.dat";
    double *data = readfile(fn);

    double llh = logLikelihood_unbinned(data, pdf);
    double llh1 = logLikelihood_unbinned(data, pdf, 0.01, 0.01);

    cout << "llh : " << llh << endl;
    cout << "llh1 : " << llh1 << endl;

    pdf->Draw();

    return;
}

double compute_likelihood_1(double sin2theta14=0, double Dm2_41=0)
{
    vInterpolator* interp0 = getPDF0();
    TF1* pdf0 = interp0->GetTF1();

    vInterpolator* interp1 = getPDF1(sin2theta14, Dm2_41);
    TF1* pdf1 = interp1->GetTF1();

    pdf0->SetNormalized(1);
    pdf1->SetNormalized(1);

    TString fn = "Cr51_L,LoE.dat";
    double *data = readfile(fn);

    double llh = logLikelihood_unbinned(data, pdf0);
    double llh1 = logLikelihood_unbinned(data, pdf1);

    cout << "llh : " << llh << endl;
    cout << "llh1 : " << llh1 << endl;

    /*
    pdf0->SetNpx(10000);
    pdf1->SetNpx(10000);
    cout << "integral pdf0 = " << pdf0->Integral(2, 19) << endl;
    cout << "integral pdf1 = " << pdf1->Integral(2, 19) << endl;
    pdf0->Draw();
    pdf1->Draw("same");
    */

    return llh1;
}

double compute_likelihood(double* data, double sin2theta14=0, double Dm2_41=0)
{
    vInterpolator* interp1 = getPDF1(sin2theta14, Dm2_41);
    TF1* pdf1 = interp1->GetTF1();

    pdf1->SetNormalized(1);

    double llh1 = logLikelihood_unbinned(data, pdf1);

    //cerr << "llh1 : " << llh1 << endl;

    /*
    pdf0->SetNpx(10000);
    pdf1->SetNpx(10000);
    cout << "integral pdf0 = " << pdf0->Integral(2, 19) << endl;
    cout << "integral pdf1 = " << pdf1->Integral(2, 19) << endl;
    pdf0->Draw();
    pdf1->Draw("same");
    */

    return llh1;
}

void compute_likelihood_grid()
{
    TString fn = "Cr51_L,LoE.dat";
    double *data = readfile(fn);

    double llh = 0; 
    int Nx = 10;
    int Ny = 10;
    double xmin = 1e-3;
    double xmax = 1;
    double ymin = 1e-2;
    double ymax = 10;

    cout << "#sin2theta14" << '\t';
    cout << "Dm2_41" << '\t';
    cout << "llh" << endl;
    double sin2theta14, Dm2_41;

    for (int i=0; i<Nx; i++) {
        sin2theta14 = TMath::Power(10, 
                        (TMath::Log10(xmax) - TMath::Log10(xmin))/Nx*i + TMath::Log10(xmin)
                      );
        for (int j=0; j<Ny; j++) {
            Dm2_41 = TMath::Power(10, 
                        (TMath::Log10(ymax) - TMath::Log10(ymin))/Ny*j + TMath::Log10(ymin)
                     );

            llh = compute_likelihood(data, sin2theta14, Dm2_41);
            cout << sin2theta14 << '\t';
            cout << Dm2_41 << '\t';
            cout << llh << endl;
        }
    }
    
}

int likelihood()
{
    compute_likelihood_grid();
    //operation_check();
    return 0;
}

