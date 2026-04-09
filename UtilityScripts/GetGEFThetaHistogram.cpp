#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

using namespace std;



void FormatHistTheta(TH1D* h, int col, int lstyle)
{
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(2);
    h->SetLineStyle(lstyle);
    h->SetStats(0);
}

void GetGEFThetaHistogram(const TString& inDir = "./GEF_tree", int Emin = 5, int Emax = 26)
{
    // ---------------- Style ----------------
    gStyle->SetPalette(57);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132,"xyz");
    gStyle->SetTitleFont(132,"a");
    gStyle->SetLabelFont(132,"xyz");
    gStyle->SetOptStat(0);
    gStyle->SetGridColor(kGray+2);
    gStyle->SetTitleFontSize(0.07);
    gROOT->ForceStyle();

    if (Emax < Emin) {
        cout << "Error: Emax < Emin" << endl;
        return;
    }

    const int nE = Emax - Emin + 1;

    // ---------------- Colors ----------------
    int numColors = nE + 1;
    int color[numColors];
    
    double startRed = 1.0;
    double startGreen = 0.0; // Adjusted for orange
    double startBlue = 1.0;
    double endRed = 1.0;
    double endGreen = 0.647;
    double endBlue = 0.0;

    // Create TGraphErrors with custom colors
    for (int i = 0; i < numColors; i++) {
        double red = startRed + (i * (endRed - startRed)) / (numColors - 1);
        double green = startGreen + (i * (endGreen - startGreen)) / (numColors - 1);
        double blue = startBlue + (i * (endBlue - startBlue)) / (numColors - 1);

        // Set the line color to the custom RGB values
        color[i] = TColor::GetColor((Float_t)red,(Float_t) green,(Float_t) blue);
    }

    if (nE > 150) {
        cout << "Warning: more than 150 excitation energies, colors will repeat." << endl;
    }

    // ---------------- Histogram arrays ----------------
    TH1D* hThetaLight[nE];
    TH1D* hThetaHeavy[nE];

    for (int i = 0; i < nE; i++) {
        hThetaLight[i] = 0;
        hThetaHeavy[i] = 0;
    }

    double maxY = 0.0;
    double maxTheta = 180.0;

    // ---------------- Read one file per E* ----------------
    for (int i = 0; i < nE; i++) {

        double Estar = Emin + i;

        TString filename = Form("%s/GEFResults_Z92_A239_E%.1f_factor_1.root",
                                inDir.Data(), Estar);

        TChain* tree = new TChain("GEFtree");
        tree->Add(filename);

        if (tree->GetEntries() == 0) {
            cout << "No entries found in: " << filename << endl;
            delete tree;
            continue;
        }

        Float_t cos_theta_light = -999.0;
        Float_t cos_theta_heavy = -999.0;

        tree->SetBranchAddress("cos_theta_light", &cos_theta_light);
        tree->SetBranchAddress("cos_theta_heavy", &cos_theta_heavy);

        hThetaLight[i] = new TH1D(
            Form("hThetaLight_E%.1f", Estar),
            "^{239}U GEF fragment #theta distributions;#theta (deg);Counts",
            300, 0.0, maxTheta);

        hThetaHeavy[i] = new TH1D(
            Form("hThetaHeavy_E%.1f", Estar),
            "^{239}U GEF fragment #theta distributions;#theta (deg);Counts",
            300, 0.0, maxTheta);

        int col = color[i % 150];
        FormatHistTheta(hThetaLight[i], col, 1); // solid = light
        FormatHistTheta(hThetaHeavy[i], col, 2); // dashed = heavy

        Long64_t nEntries = tree->GetEntries();

        for (Long64_t j = 0; j < nEntries; j++) {
            tree->GetEntry(j);

            if (cos_theta_light >= -1.0 && cos_theta_light <= 1.0) {
                double thetaL = TMath::ACos(cos_theta_light) * TMath::RadToDeg();
                if (thetaL >= 0.0 && thetaL <= maxTheta)
                    hThetaLight[i]->Fill(thetaL);
            }

            if (cos_theta_heavy >= -1.0 && cos_theta_heavy <= 1.0) {
                double thetaH = TMath::ACos(cos_theta_heavy) * TMath::RadToDeg();
                if (thetaH >= 0.0 && thetaH <= maxTheta)
                    hThetaHeavy[i]->Fill(thetaH);
            }
        }

        if (hThetaLight[i]->GetMaximum() > maxY)
            maxY = hThetaLight[i]->GetMaximum();

        if (hThetaHeavy[i]->GetMaximum() > maxY)
            maxY = hThetaHeavy[i]->GetMaximum();

        delete tree;
    }

    // ---------------- Check at least one histogram exists ----------------
    bool hasHist = false;
    for (int i = 0; i < nE; i++) {
        if (hThetaLight[i] || hThetaHeavy[i]) {
            hasHist = true;
            break;
        }
    }

    if (!hasHist) {
        cout << "No histograms were created." << endl;
        return;
    }

    // ---------------- Canvas ----------------
    TCanvas* cTheta = new TCanvas("cTheta", "GEF Theta distributions", 1920, 1080);
    cTheta->cd();
    gPad->SetGrid();

    TH1D* hFrame = new TH1D(
        "hFrame",
        "^{239}U GEF fragment #theta distributions;#theta (deg);Counts",
        300, 0.0, maxTheta);

    hFrame->SetStats(0);
    hFrame->GetXaxis()->CenterTitle(true);
    hFrame->GetYaxis()->CenterTitle(true);
    hFrame->GetXaxis()->SetTitleOffset(0.85);
    hFrame->GetYaxis()->SetTitleOffset(0.85);
    hFrame->GetXaxis()->SetTitleSize(0.05);
    hFrame->GetYaxis()->SetTitleSize(0.05);
    hFrame->GetXaxis()->SetLabelSize(0.045);
    hFrame->GetYaxis()->SetLabelSize(0.045);
    hFrame->SetMinimum(0.0);
    hFrame->SetMaximum(1.10 * maxY);
    hFrame->Draw();

    for (int i = 0; i < nE; i++) {
        if (hThetaLight[i]) hThetaLight[i]->Draw("HIST SAME");
        if (hThetaHeavy[i]) hThetaHeavy[i]->Draw("HIST SAME");
    }

    // ---------------- Legends ----------------
    TH1D* firstLight = 0;
    TH1D* firstHeavy = 0;

    for (int i = 0; i < nE; i++) {
        if (!firstLight && hThetaLight[i]) firstLight = hThetaLight[i];
        if (!firstHeavy && hThetaHeavy[i]) firstHeavy = hThetaHeavy[i];
    }

    TLegend* legType = new TLegend(0.13, 0.78, 0.27, 0.88);
    legType->SetTextFont(132);
    legType->SetTextSize(0.035);
    if (firstLight) legType->AddEntry(firstLight, "Light fragment", "l");
    if (firstHeavy) legType->AddEntry(firstHeavy, "Heavy fragment", "l");
    legType->Draw();

    TLegend* legE = new TLegend(0.72, 0.15, 0.88, 0.88);
    legE->SetTextFont(132);
    legE->SetTextSize(0.025);
    legE->SetMargin(0.35);

    for (int i = 0; i < nE; i++) {
        if (hThetaLight[i]) {
            double Estar = Emin + i;
            legE->AddEntry(hThetaLight[i], Form("E* = %.1f MeV", Estar), "l");
        }
    }
    legE->Draw();

    // ---------------- Save ----------------
    /*TFile* fout = new TFile("GEF_ThetaDistributions_U239.root", "RECREATE");
    cTheta->Write();
    hFrame->Write();

    for (int i = 0; i < nE; i++) {
        if (hThetaLight[i]) hThetaLight[i]->Write();
        if (hThetaHeavy[i]) hThetaHeavy[i]->Write();
    }

    fout->Write();
    fout->Close();

    cTheta->SaveAs("GEF_ThetaDistributions_U239.pdf");*/
}