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
    f->SetMarkerSize(0.9);
    f->SetMarkerStyle(20);
    f->SetMarkerColor(kBlue);
    f->Draw(draw_opt);

    if (mf != NULL) {
        for (int i = 0; i < mf->size(); i++) {
            (mf->at(i))->Draw("P same");
        }
    }

    TLatex *latex = new TLatex(3, 175, Form("B(I) = %.1f I + %.1f", (mf->at(0))->GetParameter("p0"), (mf->at(0))->GetParameter("p1")));
    latex->SetTextSize(0.05);
    latex->Draw("same");
    c->SaveAs(filename + Form(".png"));
    c->SaveAs(filename + Form(".root"));

    // c->SaveAs(filename + Form(".C"));

    c->Close();
}

template <typename T, typename U>
void plot_TGraph_multifunction_multipoint(U *f, U *f1, U *f2, vector<T *> *mf, TString title = "", TString draw_opt = "", TString filename = "") { // title must have same pattern of TH1 constructor => "name; axisx; axisy"
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

    f->GetYaxis()->SetLimits(0, 100);
    f->GetXaxis()->SetLimits(0.0003, 0.002);
    f->SetMinimum(0);
    f->SetMinimum(100);

    // f->SetMarkerColor(4);
    f->SetMarkerSize(0.9);
    f->SetMarkerStyle(20);
    f->SetMarkerColor(kBlue);
    f->Draw(draw_opt);

    f1->SetMarkerSize(0.9);
    f1->SetMarkerStyle(20);
    f1->SetMarkerColor(kRed);
    f1->Draw("P SAME");

    f2->SetMarkerSize(0.9);
    f2->SetMarkerStyle(20);
    f2->SetMarkerColor(kViolet);
    f2->Draw("P SAME");

    f->GetYaxis()->SetLimits(0, 100);
    f->GetXaxis()->SetLimits(0.0003, 0.002);

    if (mf != NULL) {
        for (int i = 0; i < mf->size(); i++) {
            (mf->at(i))->Draw("P same");
        }
    }

    TLatex *latex = new TLatex(3, 175, Form("B(I) = %.1f I + %.1f", (mf->at(0))->GetParameter("p0"), (mf->at(0))->GetParameter("p1")));
    latex->SetTextSize(0.05);
    latex->Draw("same");
    c->Update();
    c->SaveAs(filename + Form(".png"));
    c->SaveAs(filename + Form(".root"));

    // c->SaveAs(filename + Form(".C"));

    c->Close();
}

double campo_helmholtz(double I) {
    double N = 320, r = 12.5e-2 / 2;
    double mu0 = 1.256e-6;
    return mu0 * pow(4. / 5, 3. / 2) * N * I / r;
}

vector<TGraph *> make_fit(TString file_name = "data/dados_calibracao_Hall.txt") {
    // reading a dat file with name in file_name variable
    ifstream inp;
    inp.open(file_name);
    double _x, _y, xmax, xmin;
    vector<double> x, y;

    std::string dummy1; // just to skip the line. look for other solution later
    getline(inp, dummy1);

    while (inp >> _y >> _x) {
        if (_x > xmax) xmax = _x;
        if (_x < xmin) xmin = _x;

        x.push_back(_x);
        y.push_back(_y);
    }

    for (int i = 0; i < x.size(); i++) {
        x.at(i) = campo_helmholtz(x.at(i));
        y.at(i) = y.at(i);
    }
    double offset = (xmax - xmin) * 0.1;

    TF1 *linfac = new TF1("linfac", "[0]*x+[1]", xmin - offset, xmax + offset);
    linfac->SetNpx(10000);
    TGraph *graph1 = new TGraph(x.size(), &x.at(0), &y.at(0));
    TGraph *graph2 = new TGraph(x.size(), &x.at(0), &y.at(0));

    graph1->Fit("linfac", "R");

    vector<TF1 *> functions;
    functions.push_back(linfac);

    // getting names
    Ssiz_t from = 0;
    TString path, filename;
    file_name.Tokenize(path, from, "/");
    file_name.Tokenize(filename, from, "/");

    plot_TGraph_multifunction(graph1, &functions, ";B (T);#nu (MHz)", "AP", Form("plots/") + filename);

    vector<TGraph *> graphs_vector;
    graphs_vector.push_back(graph1);
    graphs_vector.push_back(graph2);
    return graphs_vector;
}

void analise() {
    vector<TGraph *> graphs_vector;
    graphs_vector = make_fit("data/DPPH_1.txt");
    TGraph *graph1 = graphs_vector.at(1);

    double magBohr = 9.27400899e-24; // J·T-1
    double h = 6.62607015e-34;       // m2 kg / s
    double g = h * graphs_vector.at(0)->GetFunction("linfac")->GetParameter(0) * 1.e6 / magBohr;
    double g_er = h * graphs_vector.at(0)->GetFunction("linfac")->GetParError(0) * 1.e6 / magBohr;
    cout << Form("g= %g \\pm %g", g, g_er) << endl;

    vector<TF1 *> functions;
    functions.push_back(graphs_vector.at(0)->GetFunction("linfac"));

    graphs_vector = make_fit("data/DPPH_2.txt");
    TGraph *graph2 = graphs_vector.at(1);

    g = h * graphs_vector.at(0)->GetFunction("linfac")->GetParameter(0) * 1.e6 / magBohr;
    g_er = h * graphs_vector.at(0)->GetFunction("linfac")->GetParError(0) * 1.e6 / magBohr;
    cout << Form("g= %g \\pm %g", g, g_er) << endl;

    graphs_vector = make_fit("data/DPPH_3.txt");
    TGraph *graph3 = graphs_vector.at(1);

    g = h * graphs_vector.at(0)->GetFunction("linfac")->GetParameter(0) * 1.e6 / magBohr;
    g_er = h * graphs_vector.at(0)->GetFunction("linfac")->GetParError(0) * 1.e6 / magBohr;
    cout << Form("g= %g \\pm %g", g, g_er) << endl;

    plot_TGraph_multifunction_multipoint(graph1, graph2, graph3, &functions, ";B (T);#nu  (MHz)", "AP", Form("plots/colored_DPPH_senoidal"));

    graphs_vector = make_fit("data/DPPH_triangular.txt");
    g = h * graphs_vector.at(0)->GetFunction("linfac")->GetParameter(0) * 1.e6 / magBohr;
    g_er = h * graphs_vector.at(0)->GetFunction("linfac")->GetParError(0) * 1.e6 / magBohr;
    cout << Form("g= %g \\pm %g", g, g_er) << endl;
}