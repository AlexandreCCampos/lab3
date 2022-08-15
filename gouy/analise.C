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
void plot_TGraph_multifunction_mg(U *f, vector<T *> *mf, TString title = "", TString draw_opt = "", TString filename = "") { // title must have same pattern of TH1 constructor => "name; axisx; axisy"
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
    // f->SetMarkerSize(0.9);
    // f->SetMarkerStyle(20);
    // f->SetMarkerColor(kBlue);
    // f->Draw(draw_opt);

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

TF1 *make_fit(TString file_name = "data/dados_calibracao_Hall.txt") {
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
        y.at(i) = y.at(i);
    }

    TF1 *linfac = new TF1("linfac", "[0]*x+[1]", 0, 5.1);
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

    plot_TGraph_multifunction(graph1, &functions, ";I (A); B (mT)", "AP", Form("plots/") + filename);
    return linfac;
}

TGraph *suscept_calc(TString file_name = "data/amostra_2.txt", TF1 *calibration = NULL, double l = 1, double m = 1) {
    if (calibration != NULL) {
        // reading a dat file with name in file_name variable
        ifstream inp;
        inp.open(file_name);
        double _x, _y, xmax, xmin;
        vector<double> x, y;

        std::string dummy1; // just to skip the line. look for other solution later
        getline(inp, dummy1);
        getline(inp, dummy1);

        while (inp >> _x >> _y) {
            if (_x > xmax) xmax = _x;
            if (_x < xmin) xmin = _x;

            x.push_back(_x);
            y.push_back(_y);
        }

        for (int i = 0; i < x.size(); i++) {
            x.at(i) = calibration->Eval(x.at(i) * 1.e-3 / 0.1);
            x.at(i) = x.at(i) * x.at(i);
            y.at(i) = y.at(i) * 1e3;
        }
        double offset = (xmax - xmin) * 0.1;
        TF1 *linfac = new TF1("linfac", "[0]*x+[1]", 0 - offset, pow(calibration->Eval(xmax + offset), 2));
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

        plot_TGraph_multifunction(graph1, &functions, "; B^{2} (mT)^{2};m (mg)", "AP", Form("plots/") + filename);

        double chi = 2e-3 * l * linfac->GetParameter(0) / m;
        double chi_err = 2e-3 * l * linfac->GetParError(0) / m;
        cout << Form("susceptibilidade chi = %g \\pm %g", chi, chi_err) << endl;

        return graph1;

    } else {
        cout << "no calibration" << endl;
    }
}

void analise() {
    TF1 *calibration = make_fit("data/calibration.txt");
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);

    TMultiGraph *mg = new TMultiGraph();
    TGraph *graph1 = suscept_calc("data/amostra_2.txt", calibration, 9.5, 1.48);
    TGraph *graph2 = suscept_calc("data/amostra_3.txt", calibration, 11.5, 0.4840);
    TGraph *graph3 = suscept_calc("data/amostra_4.txt", calibration, 11.3, 0.4377);
    TGraph *graph4 = suscept_calc("data/amostra_5.txt", calibration, 8.4, 2.13);

    graph1->GetFunction("linfac")->SetLineColor(kRed + 3);
    graph2->GetFunction("linfac")->SetLineColor(kBlue);
    graph3->GetFunction("linfac")->SetLineColor(kViolet);
    graph4->GetFunction("linfac")->SetLineColor(kViolet + 6);
    int i = -2;
    graph1->SetMarkerColor(kRed + 3 + i);
    graph2->SetMarkerColor(kBlue + i);
    graph3->SetMarkerColor(kViolet + i);
    graph4->SetMarkerColor(kViolet + 6 + i);

    TLegend *legend = new TLegend(0.15, 0.68, .37, .85);
    legend->AddEntry(graph1->GetFunction("linfac"), "Sal de Mohr", "l");
    legend->AddEntry(graph2->GetFunction("linfac"), "Ni CH_{2}O_{4}.2HH_{2}O - Ni#alpha", "l");
    legend->AddEntry(graph3->GetFunction("linfac"), "FeCH_{2}O_{4} 2H_{2}O - Fe#alpha", "l");
    legend->AddEntry(graph4->GetFunction("linfac"), "Co [Hg (SCN)_{4}] ", "l");

    mg->SetTitle("; B^{2} (mT)^{2};m (mg)");
    mg->Add(graph1, "AP");
    mg->Add(graph2, "AP");
    mg->Add(graph3, "AP");
    mg->Add(graph4, "AP");
    mg->Draw("a");
    legend->Draw();
    c1->SaveAs("plots/mg.png");
    // vector<TF1 *> functions;
    // functions.push_back(suscept_calc("data/amostra_2.txt", calibration, 9.5, 1.48)->GetFunction("linfac"));
    // functions.push_back(suscept_calc("data/amostra_3.txt", calibration, 11.5, 0.4840)->GetFunction("linfac"));
    // functions.push_back(suscept_calc("data/amostra_4.txt", calibration, 11.3, 0.4377)->GetFunction("linfac"));
    // functions.push_back(suscept_calc("data/amostra_5.txt", calibration, 8.4, 2.13)->GetFunction("linfac"));

    // plot_TGraph_multifunction_mg(mg, &functions, "; B^{2} (mT)^{2};m (mg)", "AP", Form("plots/mg"));
}