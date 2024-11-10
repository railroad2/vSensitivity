#include <iostream>
#include "TROOT.h"
#include "TF1.h"

#include "/home/kmlee/opt/vOscillation/vClass.hh"

//TF1* getPDF0()
vInterpolator* getPDF0()
{

    vDetectorCylinder *det = new vDetectorCylinder();
    det->SetRadius(7.5);
    det->SetHeight(15);

    double Lmin = 1.8;
    double Lmax = 20.;

    int npts = 100;

    vector<double> x, y; 

    for (int i=0; i<npts; i++) {
        double L = (Lmax - Lmin)/npts * i + Lmin;
        x.push_back(L);
        y.push_back(det->GetLDistribution(L, 10, 0, 1001, true)/ 0.747);
        cout << L << ' ' << y[i] << endl;
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

    cout << nlines << " are read." << endl;

    infile.close();

    return &(arr[0]);
}

double logLikelihood_unbinned(double* data, TF1* pdf, double theta_14=0, double Dm2_41=0)
{
    double llh = 0;
    double lh = 0;
    int nbin = 100;

    vSterile* vst = new vSterile();
    vst->Load4StdData();

    vst->Set4theta(theta_14, 0, 0); // (theta_14, theta_24, theta_34)
    vst->Set4Dm2(Dm2_41);

    for (int i=0; i<nbin; i++) {
        double tmp = pdf->Eval(data[i]/1000.);
        if (tmp <= 0) {
            tmp = 1e-31;
        }
        double osc = vst->GetProbability(data[i]/1000., 1, "e", "e");
        cout << data[i] << '\t'  <<  tmp << '\t' << osc << endl;
        llh += TMath::Log(tmp*osc);
    }
    
    return llh;
}

void likelihood()
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

int main()
{
    likelihood();
    return 0;
}

