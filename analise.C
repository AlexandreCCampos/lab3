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
    double x, y;
    int N = 0;

    inp.open(path + filename);

    std::string dummy1; // just to skip the line. look for other solution later
    getline(inp, dummy1);
    getline(inp, dummy1);
    getline(inp, dummy1);

    while (inp >> x >> y) {

        if (inp.eof()) break;

        x_array.push_back(x / 2);
        y_array.push_back(y);
        N++;
    }
    cout << "N = " << N << endl;

    inp.close();
    // achando picos
    vector<vector<double>> peaks = find_peaksn(&y_array, 200, 5);
    cout << "numero de picos = " << (peaks.at(0)).size() << endl;

    vector<double> theta_peaks;
    cout << Form("Peaks (i, x, y)") << endl;
    for (int j = 0; j < (peaks.at(0)).size(); j++) {
        cout << Form("(%d, %f, %f)", (int)peaks.at(0).at(j), x_array[(int)peaks.at(0).at(j)], peaks.at(1).at(j)) << endl;

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
    mg->SetTitle(";#theta (^{#circ});Contagens/#theta");
    // mg->GetYaxis()->SetMaxDigits(1);
    mg->GetYaxis()->SetLabelSize(0.03);

    TCanvas *c = new TCanvas("c", " ", 0, 0, 800, 600);

    mg->Add(gr);
    mg->Add(gr_peaks, "P");
    mg->Draw("A");

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

void analise() {
    read_TGraph("DRX/", "drx_Cu_LiF_abs_Ni_24kV_60s.dat");
    read_TGraph("DRX/", "drx_Cu_LiF_24kV_60s.dat");
}