#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TFile.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

void FormatHistTheta(TH1D* h, int col, int lstyle)
{
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetLineWidth(2);
    h->SetLineStyle(lstyle);
    h->SetStats(0);
}

void GetFissionFragmentsEmissionAngle(const TString& inDir = "./238U_dp_results/Event_output", int Emin = 5, int Emax = 26)
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

    double startRed   = 0.0;
    double startGreen = 0.0;
    double startBlue  = 1.0;
    double endRed     = 0.0;
    double endGreen   = 0.90;
    double endBlue    = 0.90;

    for (int i = 0; i < numColors; i++) {
        double red   = startRed   + (i * (endRed   - startRed))   / (numColors - 1);
        double green = startGreen + (i * (endGreen - startGreen)) / (numColors - 1);
        double blue  = startBlue  + (i * (endBlue  - startBlue))  / (numColors - 1);
        color[i] = TColor::GetColor((Float_t)red, (Float_t)green, (Float_t)blue);
    }

    int numColors_darker = nE + 1;
    int color_darker[numColors];

    double startRed_darker   = 0.9;
    double startGreen_darker = 0.0;
    double startBlue_darker  = 0.0;
    double endRed_darker     = 0.95;
    double endGreen_darker   = 0.0;
    double endBlue_darker    = 0.95;

    for (int i = 0; i < numColors_darker; i++) {
        double red_darker   = startRed_darker   + (i * (endRed_darker   - startRed_darker))   / (numColors_darker - 1);
        double green_darker = startGreen_darker + (i * (endGreen_darker - startGreen_darker)) / (numColors_darker - 1);
        double blue_darker  = startBlue_darker  + (i * (endBlue_darker  - startBlue_darker))  / (numColors_darker - 1);
        color_darker[i] = TColor::GetColor((Float_t)red_darker, (Float_t)green_darker, (Float_t)blue_darker);
    }

    // ---------------- Histogram arrays ----------------
    TH1D* hThetaLight[nE];
    TH1D* hThetaHeavy[nE];

    for (int i = 0; i < nE; i++) {
        hThetaLight[i] = 0;
        hThetaHeavy[i] = 0;
    }

    double maxY = 0.0;
    double maxTheta = 30.0;

    // ---------------- Read one text file per E* ----------------
    for (int i = 0; i < nE; i++) {

        double Estar = Emin + i;

        TString filename;
        if (Estar < 10.0)
            filename = Form("%s/output_event_generator_238U_dp_HRf_excEn0%.1fMeV_recoil.txt",
                            inDir.Data(), Estar);
        else
            filename = Form("%s/output_event_generator_238U_dp_HRf_excEn%.1fMeV_recoil.txt",
                            inDir.Data(), Estar);

        ifstream file(filename.Data());
        if (!file.is_open()) {
            cout << "Cannot open file: " << filename << endl;
            continue;
        }

        hThetaLight[i] = new TH1D(
            Form("hThetaLight_E%.1f", Estar),
            "^{239}U G4Beamline fragment #theta_{lab} distributions;#theta_{lab} (deg);Counts",
            150, 0.0, maxTheta);

        hThetaHeavy[i] = new TH1D(
            Form("hThetaHeavy_E%.1f", Estar),
            "^{239}U G4Beamline fragment #theta_{lab} distributions;#theta_{lab} (deg);Counts",
            150, 0.0, maxTheta);

        int col = color[i];
        int col_darker = color_darker[i];
        FormatHistTheta(hThetaLight[i], col, 1); // solid = light
        FormatHistTheta(hThetaHeavy[i], col_darker, 1); // dashed = heavy

        string line;

        // skip header
        getline(file, line);

        while (getline(file, line)) {

            if (line.empty()) continue;
            if (line[0] == '#') continue;

            stringstream ss(line);

            double x, y, z;
            double px, py, pz;
            double t;
            long long PDGid;
            int EventID, TrackID, ParentID;
            double Weight;

            ss >> x >> y >> z
               >> px >> py >> pz
               >> t >> PDGid >> EventID >> TrackID >> ParentID >> Weight;

            if (ss.fail()) continue;

            double p = sqrt(px*px + py*py + pz*pz);
            if (p <= 0.0) continue;

            double cosTheta = pz / p;
            if (cosTheta >  1.0) cosTheta =  1.0;
            if (cosTheta < -1.0) cosTheta = -1.0;

            double theta = TMath::ACos(cosTheta) * TMath::RadToDeg();

            if (theta < 0.0 || theta > maxTheta) continue;

            if (EventID % 2 == 1)
                hThetaLight[i]->Fill(theta);
            else
                hThetaHeavy[i]->Fill(theta);
        }

        file.close();

        if (hThetaLight[i]->GetMaximum() > maxY)
            maxY = hThetaLight[i]->GetMaximum();

        if (hThetaHeavy[i]->GetMaximum() > maxY)
            maxY = hThetaHeavy[i]->GetMaximum();
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
    TCanvas* cTheta = new TCanvas("cTheta", "G4Beamline Theta distributions", 1920, 1080);
    cTheta->cd();
    gPad->SetGrid();

    TH1D* hFrame = new TH1D(
        "hFrame",
        "^{239}U fragment #theta_{lab} distributions;#theta_{lab} (deg);Probability density (deg^{-1})",
        100, 0.0, 20);

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
    hFrame->SetMaximum(0.06);
    hFrame->Draw();

    for (int i = 0; i < nE; i++) {
        double integralL = hThetaLight[i]->Integral();
        double integralH = hThetaHeavy[i]->Integral();

        if (integralL > 0)
            hThetaLight[i]->Scale(1.0 / integralL);

        if (integralH > 0)
            hThetaHeavy[i]->Scale(1.0 / integralH);

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

    TLegend* legType = new TLegend(0.13, 0.78, 0.35, 0.88);
    legType->SetTextFont(132);
    legType->SetTextSize(0.035);
    if (firstLight) legType->AddEntry(firstLight, "Light fragment", "l");
    if (firstHeavy) legType->AddEntry(firstHeavy, "Heavy fragment", "l");
    legType->Draw();

    TLegend* legE = new TLegend(0.78, 0.15, 0.95, 0.88);
    legE->SetTextFont(132);
    legE->SetTextSize(0.025);
    legE->SetMargin(0.35);

    for (int i = 0; i < nE; i++) {
        if (hThetaLight[i] && hThetaHeavy[i]) {
            double Estar = Emin + i;

            legE->AddEntry(hThetaLight[i], " ", "l");
            legE->AddEntry(hThetaHeavy[i], Form("E* = %.1f MeV", Estar), "l");
        }
    }
    legE->Draw();

    cTheta->SaveAs("G4BL_ThetaDistributions_U239.pdf");

    // ---------------- Save ----------------
    /*
    TFile* fout = new TFile("G4BL_ThetaDistributions_U239.root", "RECREATE");
    cTheta->Write();
    hFrame->Write();

    for (int i = 0; i < nE; i++) {
        if (hThetaLight[i]) hThetaLight[i]->Write();
        if (hThetaHeavy[i]) hThetaHeavy[i]->Write();
    }

    fout->Write();
    fout->Close();

    
    */
}