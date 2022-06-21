#pragma once

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
#include "utility.h"
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

void plot_hist(TH1 *h, TString output_name = "") {
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "", 800, 600);
    h->GetYaxis()->SetMaxDigits(2);
    h->GetYaxis()->SetLabelSize(0.03);

    h->Draw();
    c->SaveAs(output_name + Form(".png"));
    c->Close();
}

template <typename T>
void print_vector(vector<T> *v) {
    cout << "v = [";
    for (int k = 0; k < v->size() - 1; k++) {
        cout << Form("%f,", v->at(k));
        if (k == v->size() - 2) cout << Form("%f]", v->at(k + 1)) << endl;
    }
}

// factorial of n
unsigned int fact(unsigned int n) {
    unsigned int ret = 1;
    for (unsigned int i = 1; i <= n; ++i) {
        ret *= i;
    }
    return ret;
}

// statistic from a vector - return a vector with results: mean, stddev, variance, min value and max value
std::vector<double> get_stat(std::vector<double> *sample, TString sample_name = "", bool print = 1) {

    int N = 0;
    double mean, stddev, variance, x_max = 0, x_min = 1.e6, sum1 = 0, sum2 = 0; // this values for max and min varies to each sample

    for (int i = 0; i < sample->size(); i++) {
        N++;
        x_max = max(sample->at(i), x_max);
        x_min = min(sample->at(i), x_min);
        sum1 += sample->at(i);
        sum2 += sample->at(i) * sample->at(i);
    }

    mean = sum1 / N;
    variance = sum2 / N - ((sum1 / N) * (sum1 / N));
    stddev = sqrt(variance);

    if (print) {
        // N must be equal to sample->size() and final i so
        cout << "---------------------------------------------" << endl;
        cout << "estat for " << sample_name << endl;
        cout << "N = " << N << endl;
        cout << "sample->size() = " << sample->size() << endl;
        cout << "mean = " << mean << endl;
        cout << "stddev = " << stddev << endl;
        cout << "variance = " << variance << endl;
        cout << "x_max = " << x_max << endl;
        cout << "x_min = " << x_min << endl;
        cout << "---------------------------------------------" << endl;
    }

    std::vector<double> results;
    results.push_back(mean);
    results.push_back(stddev);
    results.push_back(variance);
    results.push_back(x_min);
    results.push_back(x_max);

    return results;
}

void read_TGraph(TString path, TString output_path = "plots/") {
    ifstream inp;

    std::vector<double> x_array, y_array;
    double x;
    int N = 0;

    cout << "opening file: " << path << endl;
    inp.open(path);

    while (inp >> x) {

        x_array.push_back(N);
        y_array.push_back(x);
        N++;
        if (inp.eof()) break;
    }

    inp.close();
    TGraph *gr_peaks = new TGraph(x_array.size(), &x_array.at(0), &y_array.at(0));

    Ssiz_t from = 0;
    TString sample_name, dummy;
    path.Tokenize(dummy, from, "/");
    path.Tokenize(dummy, from, "/");
    // path.Tokenize(dummy, from, "/");
    from = 0;
    dummy.Tokenize(sample_name, from, ".");

    TCanvas *c = new TCanvas("c1", "", 800, 600);
    gr_peaks->Draw();
    c->SaveAs(output_path + sample_name + Form(".png"));
    c->Close();
}

template <typename T, typename U>
void plot_TGraph_multifunction(U *f, vector<T *> *mf, TString title = "", TString draw_opt = "", TString filename = "", double leg_pos_x = 0, double leg_pos_y = 0) { // title must have same pattern of TH1 constructor => "name; axisx; axisy"
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetGrid();
    gStyle->SetOptFit(0);
    // gStyle->SetMaxDigitsY(2);
    f->GetYaxis()->SetMaxDigits(2);
    f->GetYaxis()->SetLabelSize(0.03);
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
    f->SetMarkerSize(0.5);
    f->SetMarkerStyle(20);
    // f->SetMarkerSize(1.2);
    // f->SetMarkerStyle(4);
    // f->SetMarkerColor(kBlue);
    f->Draw(draw_opt);

    TLatex *text = new TLatex(leg_pos_x, leg_pos_y, Form("y = ax + b"));
    TLatex *text1;
    TLatex *text2;
    text->SetTextSize(0.05);

    if (mf != NULL) {
        for (int i = 0; i < mf->size(); i++) {
            (mf->at(i))->Draw("P same");
        }

        if (leg_pos_x != 0 && leg_pos_y != 0) {
            TF1 *f1 = f->GetFunction("fitlinear");
            text1 = new TLatex(leg_pos_x, leg_pos_y - (leg_pos_y * 0.1), Form("a = (%.3f #pm %.3f)", f1->GetParameter(0), f1->GetParError(0)));
            text2 = new TLatex(leg_pos_x, leg_pos_y - 2 * (leg_pos_y * 0.1), Form("b = (%.1f #pm %.1f)", f1->GetParameter(1), f1->GetParError(1)));
            text1->SetTextSize(0.05);
            text2->SetTextSize(0.05);
            text->Draw("same");
            text1->Draw("same");
            text2->Draw("same");
        }
    }

    c->SaveAs(filename + Form(".png"));
    // c->SaveAs(filename + Form(".root"));
    //  c->SaveAs(filename + Form(".C"));

    c->Close();
}

bool vector_contain(int a, vector<int> *v) {
    bool t = true;
    for (int k = 0; k < v->size(); k++) {
        if (a == v->at(k)) t = false;
    }
    // cout << "t = " << t << endl;
    return t;
}