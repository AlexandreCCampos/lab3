// find peaks in DRX spectre
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

template <typename T, typename U>
void plot_TGraph_multifunction(U *f, vector<T *> *mf, TString title = "", TString draw_opt = "", TString filename = "") { // title must have same pattern of TH1 constructor => "name; axisx; axisy"
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetGrid();

    // getting names
    Ssiz_t from = 0;
    TString gtitle, yaxistitle, xaxistitle;
    title.Tokenize(gtitle, from, ";");
    title.Tokenize(xaxistitle, from, ";");
    title.Tokenize(yaxistitle, from, ";");
    f->SetTitle(gtitle);
    f->GetXaxis()->SetTitle(xaxistitle);
    f->GetYaxis()->SetTitle(yaxistitle);
    // f->SetMarkerColor(4);
    f->SetMarkerSize(1);
    f->SetMarkerStyle(20);
    f->Draw(draw_opt);

    if (mf != NULL) {
        for (int i = 0; i < mf->size(); i++) {
            (mf->at(i))->Draw("P same");
        }
    }

    TLatex *latex = new TLatex(1.5, 110, Form("B(I) = %.1f I + %.1f", (mf->at(0))->GetParameter("p0"), (mf->at(0))->GetParameter("p1")));
    latex->SetTextSize(0.05);
    latex->Draw("same");
    c->SaveAs(filename + Form(".png"));
    c->SaveAs(filename + Form(".root"));

    // c->SaveAs(filename + Form(".C"));

    c->Close();
}

void make_fit_magnet_calibration(TString file_name = "data/dados_calibracao_Hall.txt", TString h_title = "") {
    // reading a dat file with name in file_name variable
    ifstream inp;
    inp.open(file_name);
    double _x, _y, xmax, xmin;
    vector<double> x, y;

    std::string dummy1; // just to skip the line. look for other solution later
    getline(inp, dummy1);

    while (inp >> _x >> _y) {
        if (_x > xmax) xmax = _x;
        if (_x < xmin) xmin = _x;

        x.push_back(_x);
        y.push_back(_y);
    }

    for (int i = 0; i < x.size(); i++) {
        x.at(i) = x.at(i) * 1.e-3 / 0.1;
        y.at(i) = y.at(i) * (-1);
    }

    TF1 *linfac = new TF1("linfac", "[0]*x+[1]", 0, 1.25);
    linfac->SetNpx(10000);
    TGraph *graph1 = new TGraph(x.size(), &x.at(0), &y.at(0));

    graph1->Fit("linfac", "R");

    vector<TF1 *> functions;
    functions.push_back(linfac);

    // getting names
    Ssiz_t from = 0;
    TString path, filename;
    file_name.Tokenize(path, from, "/");
    file_name.Tokenize(filename, from, "/");

    plot_TGraph_multifunction(graph1, &functions, h_title, "AP*", Form("plots/") + filename);
}

void make_fit(TString file_name = "data/dados_calibracao_Hall.txt", TString h_title = "") {
    // reading a dat file with name in file_name variable
    ifstream inp;
    inp.open(file_name);
    double _x, _y, xmax, xmin;
    vector<double> x, y;

    std::string dummy1; // just to skip the line. look for other solution later
    getline(inp, dummy1);

    while (inp >> _y >> _x) {
        _x = _x * 10;
        if (_x > xmax) xmax = _x;
        if (_x < xmin) xmin = _x;

        x.push_back(_x);
        y.push_back(_y);
    }

    for (int i = 0; i < x.size(); i++) {
        x.at(i) = x.at(i);
        y.at(i) = y.at(i);
    }

    double offset = (xmax - xmin) * 0.1;

    TF1 *linfac = new TF1("linfac", "[0]*x+[1]", xmin - offset, xmax + offset);
    linfac->SetNpx(10000);
    TGraph *graph1 = new TGraph(x.size(), &x.at(0), &y.at(0));

    graph1->Fit("linfac", "R");

    vector<TF1 *> functions;
    functions.push_back(linfac);

    // getting names
    Ssiz_t from = 0;
    TString path, filename;
    file_name.Tokenize(path, from, "/");
    file_name.Tokenize(filename, from, "/");

    plot_TGraph_multifunction(graph1, &functions, h_title, "AP*", Form("plots/") + filename);
}

void analise() {
    // make_fit_magnet_calibration("data/dados_calibracao_Hall.txt", ";I (A); B (mT)");
    make_fit("data/zinc_0.25A.txt", ";I (mA); U_{H} (mV)");
    make_fit("data/zinc_0.5A.txt", ";I (mA); U_{H} (mV)");
    make_fit("data/zinc_0.75A.txt", ";I (mA); U_{H} (mV)");
    make_fit("data/zinc_1A.txt", ";I (mA); U_{H} (mV)");
}