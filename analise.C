#include "TArrow.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <fstream>
#include <iostream>

using namespace std;

template <typename T>
vector<vector<T>> find_peaksn(vector<T> *x, double height, int order = 1) {
    int peak_found = 0;
    vector<vector<double>> peaks;
    vector<double> peak_x, peak_y;
    bool der_inf = true, der_sup = true;
    for (int i = order; i < (x->size() - order); i++) {
        for (int j = 1; j < (order + 1); j++) {
            der_inf = der_inf && (x->at(i) > x->at(i - j));
            der_sup = der_sup && (x->at(i) > x->at(i + j));
        }
        if (der_inf && der_sup && x->at(i) >= height) {
            peak_found++;
            peak_x.push_back(i);
            peak_y.push_back(x->at(i));
        }
        der_inf = true;
        der_sup = true;
    }
    peaks.push_back(peak_x);
    peaks.push_back(peak_y);
    return peaks;
}

TGraph read_TGraph(TString path, TString filename) {
    ifstream inp;

    std::vector<double> x_array, y_array;
    double x, y, x_bremsstrahlung_min = 10e3, y_bremsstrahlung_min = 10e3;
    int N = 0;

    inp.open(path + filename);

    std::string dummy1; // just to skip the line. look for other solution later
    getline(inp, dummy1);
    getline(inp, dummy1);
    getline(inp, dummy1);

    while (inp >> x >> y) {

        if (inp.eof()) break;

        if (x / 2 > 0 && x / 2 < 10 && y_bremsstrahlung_min > y) {
            x_bremsstrahlung_min = x / 2;
            y_bremsstrahlung_min = y;
        }
        x_array.push_back(x / 2);
        y_array.push_back(y);
        N++;
    }
    cout << "N = " << N << endl;
    cout << "\\theta_brehmsstralum_min = " << x_bremsstrahlung_min << endl;

    // change it later
    double dLiF = 201e-12, hc = 1.23984198e-6, q = 1.6e-19;

    double angle_rad = (x_bremsstrahlung_min)*M_PI / 180;
    double B_min_energy = hc / (2 * dLiF * sin(angle_rad) * 1.e3);
    double B_min_energy_error = hc * cos(angle_rad) * 0.5 * M_PI / 180 / (2 * dLiF * sin(angle_rad) * sin(angle_rad) * 1.e3);
    cout << "\\theta_bremsstrahlung_min_energy = " << B_min_energy << ("+/-") << B_min_energy_error << (" kV") << endl;

    inp.close();
    // achando picos
    vector<vector<double>> peaks = find_peaksn(&y_array, 200, 5);
    cout << "numero de picos = " << (peaks.at(0)).size() << endl;

    vector<double> theta_peaks;
    cout << Form("Peaks (i, x, y, E, dE)") << endl;
    for (int j = 0; j < (peaks.at(0)).size(); j++) {
        cout << Form("(%d, %f, %f, %f, %f)", (int)peaks.at(0).at(j), x_array[(int)peaks.at(0).at(j)], peaks.at(1).at(j), hc / (2 * dLiF * sin(x_array[(int)peaks.at(0).at(j)] * M_PI / 180) * 1.e3), hc * cos(x_array[(int)peaks.at(0).at(j)] * M_PI / 180) * 0.5 * M_PI / 180 / (2 * dLiF * sin(x_array[(int)peaks.at(0).at(j)] * M_PI / 180) * sin(x_array[(int)peaks.at(0).at(j)] * M_PI / 180) * 1.e3)) << endl;

        theta_peaks.push_back(x_array[(int)peaks.at(0).at(j)]);
    }

    int nPeaks = (peaks.at(0)).size();
    // fim achando picos

    TGraph *gr = new TGraph(x_array.size(), &x_array.at(0), &y_array.at(0));
    TGraph *gr_peaks = new TGraph((peaks.at(0)).size(), &theta_peaks.at(0), &peaks.at(1).at(0));

    gr_peaks->SetMarkerColor(kBlue);
    gr_peaks->SetMarkerStyle(kCircle);

    auto mg = new TMultiGraph();
    mg->SetName("mg");
    mg->SetTitle(";2#theta (^{#circ});Contagens/2#theta");

    TCanvas *c = new TCanvas("c", " ", 0, 0, 800, 600);

    mg->Add(gr);
    mg->Add(gr_peaks, "P");
    mg->Draw("A");
    mg->GetYaxis()->SetMaxDigits(2);
    mg->GetYaxis()->SetLabelSize(0.03);

    c->SaveAs(Form("plots/") + filename + Form(".png"));
    c->Close();
    // finding peak: delete later
    // print_vector(&peaks.at(0));
    // print_vector(&peaks.at(1));
    // TF1 *gauss1 = new TF1("gauss1", "gaus(0)", peaks.at(0).at((peaks.at(0)).size() - 1) - 75, peaks.at(0).at((peaks.at(0)).size() - 1) + 75);
    // TF1 *gauss2 = new TF1("gauss2", "gaus(0)", peaks.at(0).at((peaks.at(0)).size() - 2) - 75, peaks.at(0).at((peaks.at(0)).size() - 2) + 75);

    // gauss1->SetParameters(1, peaks.at(0).at((peaks.at(0)).size() - 1), 50);
    // gauss2->SetParameters(1, peaks.at(0).at((peaks.at(0)).size() - 2), 50);

    // gr_peaks->Fit(gauss1, "R");
    // gr_peaks->Fit(gauss2, "R");
    // double Delta_E = gauss1->GetParameter("Mean") - gauss2->GetParameter("Mean");
    // double mean_FWHM = sqrt(2 * log(2)) * (gauss1->GetParameter("Sigma") + gauss2->GetParameter("Sigma"));

    // cout << "separation = " << Delta_E / mean_FWHM << endl;

    //

    return *gr_peaks;
}

void superpose(TString path, TString filename1, TString filename2) {

    TGraph *gr1 = new TGraph(path + filename1, "%lg %lg", "");
    gr1->SetLineColor(kRed);
    TGraph *gr2 = new TGraph(path + filename2, "%lg %lg", "");
    gr2->SetLineColor(kBlue);

    TArrow *kb = new TArrow(41.8, 0.008e6, 41.8, 0.006e6, 0.01, "|>");
    TLatex *tkb = new TLatex(40.3, 0.0082e6, "K_{#beta}");
    tkb->SetTextSize(0.03);
    // tkb->SetTextColor(kBlue);

    kb->SetAngle(40);
    kb->SetLineWidth(2);
    // kb->SetLineColor(kBlue);
    // kb->SetFillColor(kBlue);
    TArrow *ka = new TArrow(60, 0.014e6, 50, 0.014e6, 0.01, "|>");
    TLatex *tka = new TLatex(62, 0.0138e6, "K_{#alpha}");
    tka->SetTextSize(0.03);
    // tka->SetTextColor(kBlue);

    ka->SetAngle(40);
    ka->SetLineWidth(2);
    // ka->SetLineColor(kBlue);
    // ka->SetFillColor(kBlue);
    TArrow *kbrehmmin = new TArrow(11.5, 0.0042e6, 11.5, 0.0022e6, 0.01, "|>");
    TLatex *tbrehmmin = new TLatex(10.0, 0.0044e6, "B_{min}");
    tbrehmmin->SetTextSize(0.03);
    // tbrehmmin->SetTextColor(kBlue);

    kbrehmmin->SetAngle(40);
    kbrehmmin->SetLineWidth(2);
    // kbrehmmin->SetLineColor(kBlue);
    // kbrehmmin->SetFillColor(kBlue);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1, "AL");
    mg->Add(gr2, "AL");
    mg->SetTitle(";2#theta (^{#circ});Contagens/2#theta");

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    // gr1->Draw("APL");
    // gr2->Draw("SAMEAPL");
    mg->Draw("A");
    kb->Draw();
    ka->Draw();
    tka->Draw();
    tkb->Draw();
    tbrehmmin->Draw();
    kbrehmmin->Draw();
    mg->GetYaxis()->SetMaxDigits(2);
    mg->GetYaxis()->SetLabelSize(0.03);
    c1->SaveAs("plots/overlaid_LiF_DRX.png");
    c1->Close();
}

void analise() {
    read_TGraph("DRX/", "drx_Cu_LiF_abs_Ni_24kV_60s.dat");
    read_TGraph("DRX/", "drx_Cu_LiF_24kV_60s.dat");
    superpose("DRX/", "drx_Cu_LiF_abs_Ni_24kV_60s.dat", "drx_Cu_LiF_24kV_60s.dat");
}