// print log to file in term: SomeCommand 2>&1 | tee SomeFile.txt

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
#include "TLine.h"
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

void superpose(TH1 *h1, TH1 *h2, TString filename) {

    // TArrow *kb = new TArrow(41.8, 0.008e6, 41.8, 0.006e6, 0.01, "|>");
    // TLatex *tkb = new TLatex(40.3, 0.0082e6, "K_{#beta}");
    // tkb->SetTextSize(0.03);
    // // tkb->SetTextColor(kBlue);

    // kb->SetAngle(40);
    // kb->SetLineWidth(2);
    // // kb->SetLineColor(kBlue);
    // // kb->SetFillColor(kBlue);
    // TArrow *ka = new TArrow(60, 0.014e6, 50, 0.014e6, 0.01, "|>");
    // TLatex *tka = new TLatex(62, 0.0138e6, "K_{#alpha}");
    // tka->SetTextSize(0.03);
    // // tka->SetTextColor(kBlue);

    // ka->SetAngle(40);
    // ka->SetLineWidth(2);
    // // ka->SetLineColor(kBlue);
    // // ka->SetFillColor(kBlue);
    // TArrow *kbrehmmin = new TArrow(11.5, 0.0042e6, 11.5, 0.0022e6, 0.01, "|>");
    // TLatex *tbrehmmin = new TLatex(10.0, 0.0044e6, "B_{min}");
    // tbrehmmin->SetTextSize(0.03);
    // // tbrehmmin->SetTextColor(kBlue);

    // kbrehmmin->SetAngle(40);
    // kbrehmmin->SetLineWidth(2);
    // // kbrehmmin->SetLineColor(kBlue);
    // // kbrehmmin->SetFillColor(kBlue);
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    h1->SetLineColor(kGreen);

    // h1->DrawNormalized();
    // h2->DrawNormalized("same");
    h1->SetMaximum(5000);
    h1->Draw();
    h2->Draw("same");

    c1->SaveAs("plots/overlaid_" + filename + Form(".png"));
    c1->Close();
}

TH1D *read_DRX(TString path, TString filename) {
    ifstream inp;

    // getting the sample name from path
    Ssiz_t from = 0;
    TString sample_name, dummy;
    filename.Tokenize(sample_name, from, ".");
    // end getting sample name

    std::vector<double> x_array, y_array;
    double x, y, xi, h, xf;
    int N = 0;

    inp.open(path + filename);

    std::string dummy1; // just to skip the line. look for other solution later
    getline(inp, dummy1);
    inp >> xi >> h >> xf;
    while (inp >> x) {
        x_array.push_back((xi + N * h) / 2);
        y_array.push_back(x);
        N++;

        if (inp.eof()) break;
    }
    cout << "N = " << N << endl;

    // filling histograms -make a function later
    int nbins = x_array.size();
    TH1D *h1 = new TH1D("h1", ";#theta (^{#circ});Contagens/#theta", nbins, xi / 2, xf / 2);

    for (int k = 1; k < (x_array.size() - 1); k++) {
        h1->SetBinContent(k, y_array.at(k));
    }

    // print_vector(&x_array);
    // print_vector(&y_array);

    inp.close();
    // achando picos
    vector<vector<double>> peaks = find_peaksn(&y_array, 200, 30);
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

    TCanvas *c = new TCanvas("c", " ", 0, 0, 800, 600);

    mg->Add(gr);
    mg->Add(gr_peaks, "P");
    mg->Draw("A");
    mg->GetYaxis()->SetMaxDigits(2);
    mg->GetYaxis()->SetLabelSize(0.03);

    c->SaveAs(Form("plots/") + filename + Form(".png"));
    c->Close();

    vector<double> fit_d_peaks, fit_a_peaks;
    double comp_k_alpha = 1.54184; // angstrom??????
    TF1 *gaus_peak;
    for (int k = 0; k < nPeaks; k++) {
        gaus_peak = new TF1(Form("gauss_peak_%d", k), "gaus(0)", theta_peaks.at(k) - 0.5, theta_peaks.at(k) + 0.5);
        gaus_peak->SetParameters(peaks.at(1).at(k), theta_peaks.at(k), 0.05);
        // gaus_peak->FixParameter(1, theta_peaks.at(k));

        h1->Fit(Form("gauss_peak_%d", k), "R+");
        fit_d_peaks.push_back(comp_k_alpha * 100 / (2 * sin(gaus_peak->GetParameter(1) * M_PI / 180)));
        fit_a_peaks.push_back(0);

        cout << "chi2/ndof = " << gaus_peak->GetChisquare() << "/" << gaus_peak->GetNDF() << " = " << gaus_peak->GetChisquare() / gaus_peak->GetNDF() << endl;
        double beta = (M_PI / 2) * 2.355 * gaus_peak->GetParameter(2);
        double d_cristallite = (4. / 3) * 0.9 * comp_k_alpha * 100 / (beta * cos(gaus_peak->GetParameter(1) * M_PI / 180));
        cout << Form("cristallite diameter:  d = %f pm", d_cristallite) << endl;
    }
    print_vector(&theta_peaks);
    print_vector(&fit_d_peaks);
    print_vector(&peaks.at(1));

    plot_hist(h1, Form("plots/h_") + sample_name);

    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    // c2->SetGrid();

    c2->cd();
    vector<TF1 *> net_par;
    TF1 *linear_hkl;
    TLatex *leg;
    vector<TLatex *> leg_array;

    vector<int> gambiarra_filter_miller_indices;

    // cout << "(hkl)" << endl;
    for (int h = 0; h < 5; h++) {
        for (int k = 0; k < 5; k++) {
            for (int l = 1; l < 5; l++) {
                if (vector_contain(h * h + k * k + l * l, &gambiarra_filter_miller_indices)) {
                    linear_hkl = new TF1(Form("linear_%d%d%d", h, k, l), "[0]*x", 0, 340);
                    linear_hkl->SetParameter(0, sqrt(h * h + k * k + l * l));
                    linear_hkl->SetNpx(100000);
                    linear_hkl->SetLineColor(h * h + k * k + l * l);
                    linear_hkl->SetLineWidth(1);
                    linear_hkl->SetTitle(";d (pm);a (pm)");

                    leg = new TLatex(345, linear_hkl->Eval(345), Form("(%d%d%d)", h, k, l));
                    leg->SetTextSize(0.02);
                    leg->SetTextAlign(13);
                    leg_array.push_back(leg);
                    // cout << Form("(%d%d%d)", h, k, l) << endl;
                    net_par.push_back(linear_hkl);
                    gambiarra_filter_miller_indices.push_back(h * h + k * k + l * l);
                }
            }
        }
    }

    TGraph *peaks_from_fit = new TGraph(fit_d_peaks.size(), &fit_d_peaks.at(0), &fit_a_peaks.at(0));
    peaks_from_fit->SetMarkerStyle(20);
    TLine *vertical_hkl;
    vector<TLine *> vertical_par;
    for (int k = 0; k < fit_d_peaks.size(); k++) {
        vertical_hkl = new TLine(fit_d_peaks.at(k), 0, fit_d_peaks.at(k), 1000);
        vertical_par.push_back(vertical_hkl);
    }

    TLine *horizontal_hkl;
    vector<TLine *> horizontal_par;
    cout << "------ height a (pm) ------" << endl;
    cout << Form("d \t a") << endl;
    for (int k = 0; k < net_par.size(); k++) {
        if (net_par.at(k)->Eval(fit_d_peaks.at(1)) < 1000) {
            double y = net_par.at(k)->Eval(fit_d_peaks.at(1));
            cout << Form("%f \t %f\t %f", fit_d_peaks.at(1), y, sqrt(k + 1) * 5) << endl;
            horizontal_hkl = new TLine(60, y, 340, y);
            horizontal_par.push_back(horizontal_hkl);
        }
    }
    net_par.at(net_par.size() - 1)->SetMaximum(1200);
    net_par.at(net_par.size() - 1)->Draw();
    leg_array.at(net_par.size() - 1)->Draw("SAME");
    for (int k = 0; k < net_par.size() - 1; k++) {
        net_par.at(k)->Draw("SAME");
        leg_array.at(k)->Draw("SAME");
    }
    peaks_from_fit->Draw("SAMEP");
    print_vector(&fit_d_peaks);
    for (int k = 0; k < vertical_par.size(); k++) {
        vertical_par.at(k)->Draw("SAME");
    }

    for (int k = 0; k < horizontal_par.size(); k++) {
        horizontal_par.at(k)->Draw("SAME");
    }

    c2->SaveAs(Form("plots/indices_") + filename + Form(".png"));
    c2->Close();

    // calculating 2nd and 3rd order peaks
    cout << "------ second and third orders ------" << endl;
    cout << "1st \t\t 2nd \t\t 3rd" << endl;
    for (int k = 0; k < fit_d_peaks.size(); k++) {
        double ang1, ang2, ang3;
        ang1 = asin(comp_k_alpha * 100 / (2 * fit_d_peaks.at(k))) * 180 / M_PI;
        ang2 = asin(2 * comp_k_alpha * 100 / (2 * fit_d_peaks.at(k))) * 180 / M_PI;
        ang3 = asin(3 * comp_k_alpha * 100 / (2 * fit_d_peaks.at(k))) * 180 / M_PI;

        cout << Form("%f \t %f \t %f", ang1, ang2, ang3) << endl;
    }

    // delete h1;
    return h1;
}

void analise() {

    read_DRX("data/", "KBr.dat");
    // read_DRX("data/", "KBrdet62.dat");

    TH1D *h_bruto;
    TH1D *h_moido;

    h_bruto = read_DRX("data/", "NaClbrutoPlano.dat");
    h_moido = read_DRX("data/", "NaClmoido.dat");
    superpose(h_bruto, h_moido, "NaCl");
}