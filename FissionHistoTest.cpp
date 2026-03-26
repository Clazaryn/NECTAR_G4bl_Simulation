#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <TH1D.h>
#include <TROOT.h>
#include <iostream>

using namespace std;

//____________________________________________________________
void FormatHist(TH1D* h, int col)
{
    h->SetLineColor(col);
    h->SetLineWidth(3);
    h->SetFillStyle(0);
    h->SetStats(0);
}

//____________________________________________________________
void FissionHistoTest()
{
    // --------------------------------------------------------
    // Global style, inspired by your GetGEFGraph.cpp
    // --------------------------------------------------------
    gStyle->SetPalette(57);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132,"xyz");
    gStyle->SetTitleFont(132,"a");
    gStyle->SetLabelFont(132,"xyz");
    gStyle->SetOptStat(0);
    gStyle->SetGridColor(kGray+2);
    gStyle->SetTitleFontSize(0.07);
    gROOT->ForceStyle();

    // --------------------------------------------------------
    // Custom colors close to your palette spirit
    // --------------------------------------------------------
  

    int DarkerBlueColor = TColor::GetColor((int)53,(int)39, (int)140);
    int DarkerBlueColor2 = TColor::GetColor((int)42,(int)31, (int)111);
    int DarkerBlueColor3 = TColor::GetColor((int)31,(int)23, (int)83);
    int BlueColor = TColor::GetColor((int)41,(int)95, (int)166);
    int LighterBlueColor = TColor::GetColor((int)56,(int)128, (int)223);
    int RedColor = TColor::GetColor((int)217,(int)48, (int)62);
    int DarkerRedColor = TColor::GetColor((int)129,(int)28, (int)37);
    int BrownColor = TColor::GetColor((int)72,(int)16, (int)21);
    int GreenColor = TColor::GetColor((int)10,(int)166, (int)122);
    int DarkerGreenColor = TColor::GetColor((int)5,(int)78, (int)57);
    int OrangeColor = TColor::GetColor((int)251,(int)131, (int)40);
    int DarkerOrangeColor = TColor::GetColor((int)192,(int)100, (int)31);
    int DarkerCyanColor = TColor::GetColor((int)27,(int)87, (int)107);
    int CyanColor = TColor::GetColor((int)41,(int)135, (int)166);
    int MagentaColor = TColor::GetColor((int)218,(int)48, (int)164);
    int DarkerMagentaColor = TColor::GetColor((int)129,(int)28, (int)97);
    int WaterGreenColor = TColor::GetColor((int)0,(int)224, (int)138);
    int LightGreenColor = TColor::GetColor((int)13,(int)224, (int)165);


    // --------------------------------------------------------
    // Build the chains
    // One chain for Silicon_top across all excitation energies
    // One chain for Silicon_bottom across all excitation energies
    // --------------------------------------------------------
    TChain* chTop = new TChain("Detector/Solarcells_top");
    TChain* chBot = new TChain("Detector/Silicon_bottom");
    TChain* chSide = new TChain("Detector/Silicon_side");

    chTop->Add("./238U_dp_results/Detector_output/Detectors_238U_dp_HRf_excEn*_recoil.root");
    chBot->Add("./238U_dp_results/Detector_output/Detectors_238U_dp_HRf_excEn*_recoil.root");
    chSide->Add("./238U_dp_results/Detector_output/Detectors_238U_dp_HRf_excEn*_recoil.root");

    cout << "Entries in Silicon_top chain    : " << chTop->GetEntries() << endl;
    cout << "Entries in Silicon_bottom chain : " << chBot->GetEntries() << endl;
    cout << "Entries in Silicon_side chain : " << chSide->GetEntries() << endl;

    Long64_t nTot = chTop->GetEntries() + chBot->GetEntries() + chSide->GetEntries();
    if (nTot == 0) {
        cout << "No entries found. Check:" << endl;
        cout << "  - path: ./238U_dp_results/Detector_output/" << endl;
        cout << "  - file names match Detectors_238U_dp_HRf_excEn*_recoil.root" << endl;
        cout << "  - tree names are really Detector/Silicon_top and Detector/Silicon_bottom" << endl;
        return;
    }

    // --------------------------------------------------------
    // Theta definition from momentum components
    // theta = acos(Pz / |P|)
    // in degrees
    // --------------------------------------------------------
    TString thetaExpr = "acos(Pz/sqrt(Px*Px+Py*Py+Pz*Pz))*180.0/TMath::Pi()";

    // --------------------------------------------------------
    // Histograms
    // Even EventID  -> heavy fragment
    // Odd  EventID  -> light fragment
    // --------------------------------------------------------
    const int nBins = 40;
    const double thetaMin = 0.0;
    const double thetaMax = 20.0;

    TH1D* hThetaHeavy = new TH1D("hThetaHeavy",
                                 "Heavy fragments;#theta (deg);Counts",
                                 nBins, thetaMin, thetaMax);

    TH1D* hThetaLight = new TH1D("hThetaLight",
                                 "Light fragments;#theta (deg);Counts",
                                 nBins, thetaMin, thetaMax);

    FormatHist(hThetaHeavy, RedColor);
    FormatHist(hThetaLight, BlueColor);

    // --------------------------------------------------------
    // Fill from both trees
    // even EventID = heavy
    // odd  EventID = light
    // --------------------------------------------------------
    chTop->Draw(thetaExpr + " >> hThetaHeavy", "EventID%2==0", "goff");
    chBot->Draw(thetaExpr + " >>+ hThetaHeavy", "EventID%2==0", "goff");
    chSide->Draw(thetaExpr + " >>+ hThetaHeavy", "EventID%2==0", "goff");

    chTop->Draw(thetaExpr + " >> hThetaLight", "EventID%2==1", "goff");
    chBot->Draw(thetaExpr + " >>+ hThetaLight", "EventID%2==1", "goff");
    chSide->Draw(thetaExpr + " >>+ hThetaLight", "EventID%2==1", "goff");

    cout << "Heavy entries = " << hThetaHeavy->GetEntries() << endl;
    cout << "Light entries = " << hThetaLight->GetEntries() << endl;

    // --------------------------------------------------------
    // Dummy graph for axes, like in your GetGEFGraph.cpp
    // --------------------------------------------------------
    double ymax = 1.15 * TMath::Max(hThetaHeavy->GetMaximum(), hThetaLight->GetMaximum());
    if (ymax <= 0) ymax = 1.0;

    TGraphErrors* DumbGraph = new TGraphErrors();
    DumbGraph->SetPoint(0, thetaMin, 0.0);
    DumbGraph->SetPoint(1, thetaMax, ymax);

    DumbGraph->SetTitle("Combined #theta distributions from Silicon_{top}, Silicon_{bottom} and Silicon_{side}");
    DumbGraph->SetMarkerColor(0);
    DumbGraph->SetMarkerStyle(8);
    DumbGraph->GetXaxis()->SetTitle("#theta (deg)");
    DumbGraph->GetYaxis()->SetTitle("Counts");
    DumbGraph->GetXaxis()->CenterTitle(true);
    DumbGraph->GetYaxis()->CenterTitle(true);
    DumbGraph->GetXaxis()->SetRangeUser(thetaMin, thetaMax);
    DumbGraph->GetYaxis()->SetRangeUser(0.0, ymax);
    DumbGraph->GetXaxis()->SetTitleOffset(0.8);
    DumbGraph->GetYaxis()->SetTitleOffset(0.8);
    DumbGraph->GetXaxis()->SetTitleSize(0.055);
    DumbGraph->GetYaxis()->SetTitleSize(0.055);

    // --------------------------------------------------------
    // Legend
    // --------------------------------------------------------
    TLegend* leg = new TLegend(0.64, 0.70, 0.98, 0.86);
    leg->SetTextFont(132);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->AddEntry(hThetaHeavy, "Even EventID = heavy fragment", "l");
    leg->AddEntry(hThetaLight, "Odd EventID = light fragment", "l");

    // --------------------------------------------------------
    // Main canvas with the two histograms
    // --------------------------------------------------------
    TCanvas* cTheta = new TCanvas("cTheta", "Theta distributions", 1920, 1080);
    cTheta->SetGrid();

    DumbGraph->Draw("AP");
    hThetaHeavy->Draw("HIST SAME");
    hThetaLight->Draw("HIST SAME");
    leg->Draw("SAME");

    cTheta->SaveAs("Theta_distributions_top_bottom_allExcEn.pdf");
    //cTheta->SaveAs("Theta_distributions_top_bottom_allExcEn.png");

    // --------------------------------------------------------
    // Optional: separate canvases
    // --------------------------------------------------------
    /*TCanvas* cHeavy = new TCanvas("cHeavy", "Heavy theta", 1400, 900);
    cHeavy->SetGrid();
    hThetaHeavy->Draw("HIST");
    cHeavy->SaveAs("Theta_heavy_allExcEn.pdf");

    TCanvas* cLight = new TCanvas("cLight", "Light theta", 1400, 900);
    cLight->SetGrid();
    hThetaLight->Draw("HIST");
    cLight->SaveAs("Theta_light_allExcEn.pdf");*/

    // --------------------------------------------------------
    // Save output objects in a ROOT file
    // --------------------------------------------------------
    /*TFile* fout = new TFile("Theta_distributions_top_bottom_allExcEn.root", "RECREATE");
    hThetaHeavy->Write();
    hThetaLight->Write();
    DumbGraph->Write("DumbGraph_theta");
    cTheta->Write();
    //cHeavy->Write();
    //cLight->Write();
    fout->Close();*/

    cout << "Output written to:" << endl;
    cout << "  Theta_distributions_top_bottom_allExcEn.pdf" << endl;
    //cout << "  Theta_heavy_allExcEn.pdf" << endl;
    //cout << "  Theta_light_allExcEn.pdf" << endl;
    cout << "  Theta_distributions_top_bottom_allExcEn.root" << endl;
}