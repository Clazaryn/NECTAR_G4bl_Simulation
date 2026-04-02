#include <TSystem.h>
//namespace fs = std::filesystem;

#include "TSpectrum.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TPolyLine.h"
using namespace std;  


#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include "TTree.h"
#include <string>
#include "TH2.h"
#include "TString.h"
#include <vector>
#include "TLeaf.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include <fstream>
#include <sstream>
//#include <direct.h>
//#include <bits/stdc++.h>
#include <TSystem.h>
//namespace fs = std::filesystem;
#include "TSpectrum.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <TCut.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TROOT.h>
using namespace std;  
#include <TGraph2D.h>
#include <iomanip>
#include <TApplication.h>
#include <TBrowser.h>
#include <TView.h>
#include "TView3D.h"
#include <cmath>
#include "TVirtualFFT.h"
#include <TImage.h>

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//****************************************************************************************************************************************************************************************************
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
#include <string>
#include <regex>
#include <algorithm>
#include <dirent.h>
#include <TLatex.h>
#include <TGaxis.h>
#include "det_analysis.h"

using namespace std;
//namespace fs = std::filesystem;

double* Excitationenergies;
double EexcBinSize;
int Bin_number;


//____________________________________________________________
void FormatHist(TH1D* h, int col)
{
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->SetFillStyle(1001);
    h->SetFillColor(col);
    h->SetStats(0);
}

//____________________________________________________________
// Helper find the number of files that exist and that i can loop on 
// we do that in order to not be dependent on how the previous step worked or not


bool FindFissionEnergies(const string& dirPath,
                         int& energyNumber,
                         double*& Excitationenergies)
{
    energyNumber = 0;
    Excitationenergies = nullptr;

    // =========================================================
    // First pass: count HRf files
    // =========================================================
    DIR* dir = opendir(dirPath.c_str());
    if (!dir) {
        cerr << "Error: cannot open directory " << dirPath << endl;
        return false;
    }

    int count = 0;
    struct dirent* entry;

    while ((entry = readdir(dir)) != NULL) {
        string filename = entry->d_name;

        if (filename.find("HRf") != string::npos &&
            filename.find(".root") != string::npos &&
            filename.find("excEn") != string::npos)
        {
            count++;
        }
    }

    closedir(dir);

    if (count == 0) {
        cerr << "No HRf root files found in " << dirPath << endl;
        return false;
    }

    double* tempEnergies = new double[count];

    // =========================================================
    // Second pass: extract energies
    // =========================================================
    dir = opendir(dirPath.c_str());
    if (!dir) {
        cerr << "Error: cannot reopen directory " << dirPath << endl;
        delete[] tempEnergies;
        return false;
    }

    int i = 0;
    while ((entry = readdir(dir)) != NULL) {
        string filename = entry->d_name;

        if (filename.find("HRf") == string::npos) continue;
        if (filename.find(".root") == string::npos) continue;

        size_t pos1 = filename.find("excEn");
        size_t pos2 = filename.find("MeV");

        if (pos1 == string::npos || pos2 == string::npos) continue;

        pos1 += 5; // move after "excEn"

        string energyStr = filename.substr(pos1, pos2 - pos1);

        tempEnergies[i] = atof(energyStr.c_str());
        i++;
    }

    closedir(dir);

    // =========================================================
    // Sort
    // =========================================================
    sort(tempEnergies, tempEnergies + count);

    // =========================================================
    // Remove duplicates
    // =========================================================
    int uniqueCount = 1;
    for (int j = 1; j < count; j++) {
        if (tempEnergies[j] != tempEnergies[uniqueCount - 1]) {
            tempEnergies[uniqueCount] = tempEnergies[j];
            uniqueCount++;
        }
    }

    // =========================================================
    // Final exact-size allocation
    // =========================================================
    energyNumber = uniqueCount;
    Excitationenergies = new double[energyNumber];

    for (int j = 0; j < energyNumber; j++) {
        Excitationenergies[j] = tempEnergies[j];
    }

    delete[] tempEnergies;

    return true;
}

//____________________________________________________________
// Helper get in which bin a Eex is ! 

int GetBinIndex(double Eex, double Emin, double EexcBinSize, int Bin_number)
{
    int bin = (int)((Eex - Emin) / EexcBinSize);
    if (bin < 0 || bin >= Bin_number) return -1;
    return bin;
}


//____________________________________________________________
void MakeFissionTransmissionPlot(const char* reaction_name = "238U_dp"){

    //string ReactionName = "238U_dp";
    string ReactionName = reaction_name;
    EexcBinSize = 1.0;

    // -------------------------------------------------------------------------
    // Let's find out what output file we have and loop through them 
    // -------------------------------------------------------------------------
    
    // Still to get to read the type of reaction !! 

    int energyNumber;
    Excitationenergies = nullptr;

    bool ok = FindFissionEnergies("./"+ReactionName+"_results/Det_analysis",
                                energyNumber,
                                Excitationenergies);

    energyNumber--;

    double HalfEnergyStep = (Excitationenergies[1]-Excitationenergies[0])/2.0;

    //if (!ok) return 1;

    /*for (int i = 0; i < energyNumber; i++) {
        cout << i << "  " << Excitationenergies[i] << endl;
    }*/

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
    int LighterCyanColor = TColor::GetColor((int)54,(int)176, (int)217);
    int MagentaColor = TColor::GetColor((int)218,(int)48, (int)164);
    int DarkerMagentaColor = TColor::GetColor((int)129,(int)28, (int)97);
    int WaterGreenColor = TColor::GetColor((int)0,(int)224, (int)138);
    int LightGreenColor = TColor::GetColor((int)13,(int)224, (int)165);
    int YellowColor = TColor::GetColor((int)237,(int)196,(int)59);
    int PurpleColor = TColor::GetColor((int)146,(int)10,(int)65);

    // ---------------------------------------------------------------------------------------------------------------------
    // Let's loop on the fission ".root" files and gets the number of events (and proportion) for each of the detection case  
    // ---------------------------------------------------------------------------------------------------------------------
    

    LightEjectile* ejectile = new LightEjectile();
    FissionEvent* fission = new FissionEvent();

    Int_t Z, A;
    Int_t detector_id;  // 0=primary, 1=auxillary for New detectors
    Int_t vert_strip, hor_strip;
    Double_t true_Eexc, true_Eejc, true_theta;
    Double_t recon_Eexc, recon_Eejc, recon_theta;
    Float_t meas_dE, meas_E1, meas_Eres;  // Measured Eres energy (crystal in New setup)
    Float_t meas_E2, meas_E3, meas_E4, meas_E5, meas_E6;  // Measured E2 through E6 energy (thick Si in PoP setup)
    

    Bool_t light_hit_Sidedetec;
    Bool_t light_hit_Topdetec;
    Bool_t light_hit_Bottomdetec;

    Bool_t heavy_hit_Sidedetec;
    Bool_t heavy_hit_Topdetec;
    Bool_t heavy_hit_Bottomdetec;

    Int_t light_Z, light_A;
    Int_t heavy_Z, heavy_A;

    int Full_count_side = 0;
    int Full_count_top = 0;
    int Full_count_bottom = 0;
    int Full_count_top_and_bottom = 0;
    int Full_count_top_and_side = 0;
    int Full_count_bottom_and_side = 0;

    int Full_count_no_hit = 0;
    int Full_count_all_three = 0;

    //Bin_number = (int)((Excitationenergies[energyNumber]-Excitationenergies[0])/EexcBinSize);

    Bin_number = (int)std::ceil((Excitationenergies[energyNumber] - Excitationenergies[0]) / EexcBinSize);
    double Emin = Excitationenergies[0];
    double Emax = Emin + Bin_number * EexcBinSize;

    double* entries_per_bin = new double[Bin_number];
    double* weight_per_bin  = new double[Bin_number];
    double weight = 1.0;

    for (int i = 0; i < Bin_number; i++) {
        entries_per_bin[i] = 0.0;
        weight_per_bin[i]  = 0.0;
    }



    //cout << "__from__" << Excitationenergies[0] << "__to__" << Excitationenergies[energyNumber] << "___binSize__" << EexcBinSize << "___bin_number__" 
    //     << Bin_number << "__check__" <<  Excitationenergies[energyNumber] << "__==__" << Excitationenergies[0]+(EexcBinSize*Bin_number) << endl;

    TH1D* HitRatio_Bottom = new TH1D("HitRatio_Bottom","HitRatio_Bottom",Bin_number,Emin,Emax);
    TH1D* HitRatio_Side = new TH1D("HitRatio_Side","HitRatio_Side",Bin_number,Emin,Emax);
    TH1D* HitRatio_Top = new TH1D("HitRatio_Top","HitRatio_Top",Bin_number,Emin,Emax);
    TH1D* HitRatio_TopandBottom = new TH1D("HitRatio_TopandBottom","HitRatio_TopandBottom",Bin_number,Emin,Emax);
    TH1D* HitRatio_TopandSide = new TH1D("HitRatio_TopandSide","HitRatio_TopandSide",Bin_number,Emin,Emax);
    TH1D* HitRatio_BottomandSide = new TH1D("HitRatio_BottomandSide","HitRatio_BottomandSide",Bin_number,Emin,Emax);
    

    FormatHist(HitRatio_Bottom,RedColor);
    FormatHist(HitRatio_Side,BlueColor);
    FormatHist(HitRatio_Top,GreenColor);
    FormatHist(HitRatio_TopandBottom,YellowColor);
    FormatHist(HitRatio_TopandSide,LighterCyanColor);
    FormatHist(HitRatio_BottomandSide,PurpleColor);
    
    for (int Eex_bin = 0; Eex_bin < energyNumber; Eex_bin++){

        // Let's find and open the tree associated with the reaction ! 

        TChain * tree = new TChain("events");
        if(Excitationenergies[Eex_bin]>9.999){
            tree->Add(("./"+ReactionName+"_results/Det_analysis/events_"+ReactionName+Form("_HRf_excEn%.1fMeV.root",Excitationenergies[Eex_bin])).c_str());
        }
        if(Excitationenergies[Eex_bin]<10){
            tree->Add(("./"+ReactionName+"_results/Det_analysis/events_"+ReactionName+Form("_HRf_excEn0%.1fMeV.root",Excitationenergies[Eex_bin])).c_str());
        }

        tree->SetBranchAddress("ejectile", &ejectile);
        tree->SetBranchAddress("fission", &fission);
        
        Long64_t entries = tree->GetEntries();

        // First we get the number of enetries per bin and the associated wieght of a hit event
        for (Long64_t i = 0; i < entries; i++) {
            tree->GetEntry(i);

            double Eex = ejectile->true_Eexc;
            int bin = GetBinIndex(Eex, Excitationenergies[0], EexcBinSize, Bin_number);

            if (bin >= 0) {
                entries_per_bin[bin] += 1.0;
            }
        }
        
        delete tree;
    }

    for (int i = 0; i < Bin_number; i++) {
        if (entries_per_bin[i] > 0.0) {
            weight_per_bin[i] = 100.0 / entries_per_bin[i];
        }
        else {
            weight_per_bin[i] = 0.0;
        }

        //cout << "bin " << i
        //    << "  E* in [" << Excitationenergies[0] + i * EexcBinSize
        //    << ", " << Excitationenergies[0] + (i + 1) * EexcBinSize
        //    << "]  entries = " << entries_per_bin[i]
        //    << "  weight = " << weight_per_bin[i] << endl;
    }

    for (int Eex_bin = 0; Eex_bin < energyNumber; Eex_bin++){

        // Let's find and open the tree associated with the reaction ! 

        TChain * tree = new TChain("events");
        if(Excitationenergies[Eex_bin]>9.999){
            tree->Add(("./"+ReactionName+"_results/Det_analysis/events_"+ReactionName+Form("_HRf_excEn%.1fMeV.root",Excitationenergies[Eex_bin])).c_str());
        }
        if(Excitationenergies[Eex_bin]<10){
            tree->Add(("./"+ReactionName+"_results/Det_analysis/events_"+ReactionName+Form("_HRf_excEn0%.1fMeV.root",Excitationenergies[Eex_bin])).c_str());
        }

        tree->SetBranchAddress("ejectile", &ejectile);
        tree->SetBranchAddress("fission", &fission);
        
        Long64_t entries = tree->GetEntries();


        int count_side = 0;
        int count_top = 0;
        int count_bottom = 0;
        int count_top_and_bottom = 0;
        int count_top_and_side = 0;
        int count_bottom_and_side = 0;

        int count_no_hit = 0;
        int count_all_three = 0;

        bool topHit = 0;
        bool bottomHit = 0;
        bool sideHit = 0;
        int nDetHit = 0;
        
        double Eex = 0.0;

        for (Long64_t i = 0; i < entries; i++) {

            //if (i % 5000 == 0) {
            //    double prog = 100.0 * (double)i / (double)entries;
            //    cout << "\rProgression " << prog << " %" << flush;
            //}

            tree->GetEntry(i);  

            Eex = ejectile->true_Eexc;

            int bin = GetBinIndex(Eex, Excitationenergies[0], EexcBinSize, Bin_number);
            if (bin < 0) continue;

            weight = weight_per_bin[bin];

            
            //if(ejectile->true_Eexc<Excitationenergies[Eex_bin]){
            //    cout << "__File__" <<Excitationenergies[Eex_bin] << "__entry_" << i << "__ejectile_E*__" << ejectile->true_Eexc 
            //         << "__light__hit(Top,Bottom,Side)__(" << fission->light.hit_Topdetec << ", " << fission->light.hit_Bottomdetec << ", " << fission->light.hit_Sidedetec << ") " 
            //         << "__heavy__hit(Top,Bottom,Side)__(" << fission->heavy.hit_Topdetec << ", " << fission->heavy.hit_Bottomdetec << ", " << fission->heavy.hit_Sidedetec << ") " 
            //         << endl;
            //    cout << "__File__" <<Excitationenergies[Eex_bin] << "__entry_" << i << "__ejectile_E*__" << ejectile->true_Eexc << endl;
            //}

            // Check if one fragment  or the other hit the detector
            topHit = fission->light.hit_Topdetec || fission->heavy.hit_Topdetec;
            bottomHit = fission->light.hit_Bottomdetec || fission->heavy.hit_Bottomdetec;
            sideHit = fission->light.hit_Sidedetec || fission->heavy.hit_Sidedetec;

            // Count number of detector types hit in this event
            nDetHit = (int)topHit + (int)bottomHit + (int)sideHit;

            

            if (nDetHit == 0) {
                count_no_hit++;
                Full_count_no_hit++;
            }


            else if (nDetHit == 1) {
                if (bottomHit) {
                    count_bottom++;
                    Full_count_bottom++;
                    HitRatio_Bottom->Fill(Eex, weight); // We fill with a stack up philosophy ... 
                    HitRatio_Top->Fill(Eex, weight);
                    HitRatio_Side->Fill(Eex, weight);
                    HitRatio_TopandBottom->Fill(Eex, weight);
                    HitRatio_TopandSide->Fill(Eex, weight);
                    HitRatio_BottomandSide->Fill(Eex, weight);
                }
                else if (topHit) {
                    count_top++;
                    Full_count_top++;
                    HitRatio_Top->Fill(Eex, weight);
                    HitRatio_Side->Fill(Eex, weight);
                    HitRatio_TopandBottom->Fill(Eex, weight);
                    HitRatio_TopandSide->Fill(Eex, weight);
                    HitRatio_BottomandSide->Fill(Eex, weight);
                }
                else if (sideHit) {
                    count_side++;
                    Full_count_side++;
                    HitRatio_Side->Fill(Eex, weight);
                    HitRatio_TopandBottom->Fill(Eex, weight);
                    HitRatio_TopandSide->Fill(Eex, weight);
                    HitRatio_BottomandSide->Fill(Eex, weight);
                }
            }
            else if (nDetHit == 2) {
                if (topHit && bottomHit) {
                    count_top_and_bottom++;
                    Full_count_top_and_bottom++;
                    HitRatio_TopandBottom->Fill(Eex, weight);
                    HitRatio_TopandSide->Fill(Eex, weight);
                    HitRatio_BottomandSide->Fill(Eex, weight);
                }
                else if (topHit && sideHit) {
                    count_top_and_side++;
                    Full_count_top_and_side++;
                    HitRatio_TopandSide->Fill(Eex, weight);
                    HitRatio_BottomandSide->Fill(Eex, weight);
                }
                else if (bottomHit && sideHit) {
                    count_bottom_and_side++;
                    Full_count_bottom_and_side++;
                    HitRatio_BottomandSide->Fill(Eex, weight);
                }
            }
            else if (nDetHit == 3) {
                cout << "Wait something weird happen the 3 detectors saw a hit !! " << endl;
            }

           
        }

        //HitRatio_Bottom->Fill(Excitationenergies[Eex_bin]+HalfEnergyStep,100*(count_bottom)/entries);
        //HitRatio_Top->Fill(Excitationenergies[Eex_bin]+HalfEnergyStep,100*(count_bottom+count_top)/entries);
        //HitRatio_Side->Fill(Excitationenergies[Eex_bin]+HalfEnergyStep,100*(count_bottom+count_top+count_side)/entries);
        //HitRatio_TopandBottom->Fill(Excitationenergies[Eex_bin]+HalfEnergyStep,100*(count_bottom+count_top+count_side+count_top_and_bottom)/entries);
        //HitRatio_TopandSide->Fill(Excitationenergies[Eex_bin]+HalfEnergyStep,100*(count_bottom+count_top+count_side+count_top_and_bottom+count_top_and_side)/entries);
        //HitRatio_BottomandSide->Fill(Excitationenergies[Eex_bin]+HalfEnergyStep,100*(count_bottom+count_top+count_side+count_top_and_bottom+count_top_and_side+count_bottom_and_side)/entries);

        //cout << "------------------------------" << endl;
        //cout << "___TEST___" << Excitationenergies[Eex_bin]+HalfEnergyStep << "___" << 100*(count_bottom)/entries << endl;
        //cout << "------File--E*---" << Excitationenergies[Eex_bin] << "------------" << endl;
        //cout << "\nOnly Top: " << count_top << endl;
        //cout << "Only Bottom: " << count_bottom << endl;
        //cout << "Only Side: " << count_side << endl;
        //cout << "Top + Bottom: " << count_top_and_bottom << endl;
        //cout << "Top + Side: " << count_top_and_side << endl;
        //cout << "Bottom + Side: " << count_bottom_and_side << endl;
        //cout << "All three: " << count_all_three << endl;
        //cout << "No hit: " << count_no_hit << endl;

        
        delete tree;
        tree = nullptr;   // good practice

    }




    // --------------------------------------------------------
    // Dummy graph for axes
    // --------------------------------------------------------
    double ymin = 0.0;
    double ymax = 110;
    double xmin = 0.0;
    double xmax = Excitationenergies[energyNumber]+8.0;

    TGraphErrors* DumbGraph = new TGraphErrors();
    DumbGraph->SetPoint(0, xmin, ymin);
    DumbGraph->SetPoint(1, xmax, ymax);

    DumbGraph->SetTitle("Fission's Detectors' Efficiencies");
    DumbGraph->SetMarkerColor(0);
    DumbGraph->SetMarkerStyle(8);
    DumbGraph->GetXaxis()->SetTitle(" E^{*} (MeV)");
    DumbGraph->GetYaxis()->SetTitle("Efficiency");
    DumbGraph->GetXaxis()->CenterTitle(true);
    DumbGraph->GetYaxis()->CenterTitle(true);
    DumbGraph->GetXaxis()->SetRangeUser(xmin, xmax);
    DumbGraph->GetYaxis()->SetRangeUser(ymin, ymax);
    DumbGraph->GetXaxis()->SetTitleOffset(0.8);
    DumbGraph->GetYaxis()->SetTitleOffset(0.8);
    DumbGraph->GetXaxis()->SetTitleSize(0.055);
    DumbGraph->GetYaxis()->SetTitleSize(0.055);

    // --------------------------------------------------------
    // Legend
    // --------------------------------------------------------
    TLegend* leg = new TLegend(0.7975182,0.7556818,0.9892815,0.9715909);
    leg->SetTextFont(132);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->AddEntry(HitRatio_BottomandSide, "Bottom & Side", "lf");
    leg->AddEntry(HitRatio_TopandSide, "Top & Side", "lf");
    leg->AddEntry(HitRatio_TopandBottom, "Top & Bottom", "lf");
    leg->AddEntry(HitRatio_Side, "Only Side", "lf");
    leg->AddEntry(HitRatio_Top, "Only Top", "lf");
    leg->AddEntry(HitRatio_Bottom, "Only Bottom", "lf");

    // --------------------------------------------------------
    // Main canvas with the two histograms
    // --------------------------------------------------------
    TCanvas* cFission = new TCanvas("cFission", "Fission Efficiencies", 1920, 1080);

    // Left pad: main plot
    TPad* padMain = new TPad("padMain", "padMain", 0.00, 0.00, 0.80, 1.00);
    padMain->SetLeftMargin(0.10);
    padMain->SetRightMargin(0.02);
    padMain->SetBottomMargin(0.12);
    padMain->SetTopMargin(0.08);
    padMain->Draw();

    TPad* padBar = new TPad("padBar", "padBar", 0.80, 0.00, 0.98, 0.92);
    padBar->SetLeftMargin(0.30);
    padBar->SetRightMargin(0.05);
    padBar->SetBottomMargin(0.12);
    padBar->SetTopMargin(0.08);
    padBar->Draw();

    padMain->cd();

    DumbGraph->Draw("AP");
    HitRatio_BottomandSide->Draw("HistSame");
    HitRatio_TopandSide->Draw("HistSame");
    HitRatio_TopandBottom->Draw("HistSame");
    HitRatio_Side->Draw("HistSame");
    HitRatio_Top->Draw("HistSame");
    HitRatio_Bottom->Draw("HistSame");

    padMain->SetGrid();
    gPad->RedrawAxis("G");
    gPad->RedrawAxis();

    leg->Draw("same");


    // --------------------------------------------------------
    // Now let's draw a  Side bar on the right with the FUll_count ratio
    // --------------------------------------------------------

    padBar->cd();

    double total = Full_count_bottom
             + Full_count_top
             + Full_count_side
             + Full_count_top_and_bottom
             + Full_count_top_and_side
             + Full_count_bottom_and_side
             + Full_count_no_hit;

    // optional dummy frame
    TH1F* hFrame = new TH1F("hFrame","",1,0.0,1.0);
    hFrame->SetMinimum(-1.8);
    hFrame->SetMaximum(100.0);
    hFrame->GetXaxis()->SetLabelSize(0.0);
    hFrame->GetXaxis()->SetTickLength(0.0);
    hFrame->GetXaxis()->SetAxisColor(0);
    hFrame->GetYaxis()->SetTitle("Global Efficiencies (%)");
    hFrame->GetYaxis()->CenterTitle(true);
    hFrame->GetYaxis()->SetTitleFont(132);
    hFrame->GetYaxis()->SetLabelFont(132);
    hFrame->GetYaxis()->SetTitleSize(0.12);
    hFrame->GetYaxis()->SetLabelSize(0.10);
    hFrame->GetYaxis()->SetTitleOffset(0.9);
    hFrame->Draw();

    double x1 = 0.10;
    double x2 = 0.60;

    double y0 = 0.0;

    TLatex latex;
    latex.SetTextFont(132);
    latex.SetTextAlign(12);   // left aligned
    latex.SetTextSize(0.07);

    int smallLabelIndex = 0;

    auto drawFracBox = [&](double count, int color)
    {
        double frac = 100.0 * count / total;
        TBox* b = new TBox(x1, y0, x2, y0 + frac);
        b->SetFillColor(color);
        b->SetLineColor(kBlack);
        b->Draw("same");

        // draw label to the right of the box
        double yText = y0 + 0.5 * frac;

        if (frac < 3.0) {
            yText = y0 + frac + 3.0 * smallLabelIndex;
            smallLabelIndex++;
        }
        
        latex.SetTextColor(color);
        latex.DrawLatex(x2 + 0.06, yText, Form("%.2f%%", frac));
           

        y0 += frac;
    };

    drawFracBox(Full_count_bottom, RedColor);
    drawFracBox(Full_count_top, GreenColor);
    drawFracBox(Full_count_side, BlueColor);
    drawFracBox(Full_count_top_and_bottom, YellowColor);
    drawFracBox(Full_count_top_and_side, LighterCyanColor);
    drawFracBox(Full_count_bottom_and_side, PurpleColor);

    gPad->RedrawAxis();
    // --------------------------------------------------------
    // Now we save the file 
    // --------------------------------------------------------
    cFission->SaveAs("FissionEfficiencies.pdf");
    //cFission->SaveAs("Theta_distributions_top_bottom_allExcEn.png");

    // --------------------------------------------------------
    // Create ROOT-safe canvas copy with standard colors
    // --------------------------------------------------------
TCanvas* FissionEfficiencies = new TCanvas("FissionEfficiencies",
                                           "Fission Efficiencies",
                                           1920, 1080);    FissionEfficiencies->cd();

        TH1D* HitRatio_Bottom_store = (TH1D*)HitRatio_Bottom->Clone("HitRatio_Bottom_store");
        TH1D* HitRatio_Side_store = (TH1D*)HitRatio_Side->Clone("HitRatio_Side_store");
        TH1D* HitRatio_Top_store = (TH1D*)HitRatio_Top->Clone("HitRatio_Top_store");
        TH1D* HitRatio_TopandBottom_store = (TH1D*)HitRatio_TopandBottom->Clone("HitRatio_TopandBottom_store");
        TH1D* HitRatio_TopandSide_store = (TH1D*)HitRatio_TopandSide->Clone("HitRatio_TopandSide_store");
        TH1D* HitRatio_BottomandSide_store = (TH1D*)HitRatio_BottomandSide->Clone("HitRatio_BottomandSide_store");
        
        FormatHist(HitRatio_Bottom_store, kRed+1);
        FormatHist(HitRatio_Side_store, kBlue+1);
        FormatHist(HitRatio_Top_store, kGreen+2);
        FormatHist(HitRatio_TopandBottom_store, kYellow+1);
        FormatHist(HitRatio_TopandSide_store, kCyan+1);
        FormatHist(HitRatio_BottomandSide_store, kMagenta+1);

        // Left pad
        TPad* padMain_store = new TPad("padMain_store", "padMain_store",
                                    0.00, 0.00, 0.80, 1.00);
        padMain_store->SetLeftMargin(0.10);
        padMain_store->SetRightMargin(0.02);
        padMain_store->SetBottomMargin(0.12);
        padMain_store->SetTopMargin(0.08);
        padMain_store->Draw();

        // Right pad
        TPad* padBar_store = new TPad("padBar_store", "padBar_store",
                                    0.80, 0.00, 0.98, 0.92);
        padBar_store->SetLeftMargin(0.30);
        padBar_store->SetRightMargin(0.05);
        padBar_store->SetBottomMargin(0.12);
        padBar_store->SetTopMargin(0.08);
        padBar_store->Draw();


        // ========================================================
        // MAIN PAD
        // ========================================================
        padMain_store->cd();

        TGraphErrors* DumbGraph_store = (TGraphErrors*)DumbGraph->Clone("DumbGraph_store");
        DumbGraph_store->Draw("AP");

        HitRatio_BottomandSide_store->Draw("HistSame");
        HitRatio_TopandSide_store->Draw("HistSame");
        HitRatio_TopandBottom_store->Draw("HistSame");
        HitRatio_Side_store->Draw("HistSame");
        HitRatio_Top_store->Draw("HistSame");
        HitRatio_Bottom_store->Draw("HistSame");

        padMain_store->SetGrid();
        gPad->RedrawAxis("G");
        gPad->RedrawAxis();

        // legend clone
        TLegend* leg_store = new TLegend(0.7975182,0.7556818,0.9892815,0.9715909);
        leg_store->SetTextFont(132);
        leg_store->SetTextSize(0.035);
        leg_store->SetBorderSize(1);
        leg_store->SetFillStyle(1001);
        leg_store->SetFillColor(kWhite);
        leg_store->AddEntry(HitRatio_BottomandSide_store, "Bottom & Side", "lf");
        leg_store->AddEntry(HitRatio_TopandSide_store, "Top & Side", "lf");
        leg_store->AddEntry(HitRatio_TopandBottom_store, "Top & Bottom", "lf");
        leg_store->AddEntry(HitRatio_Side_store, "Only Side", "lf");
        leg_store->AddEntry(HitRatio_Top_store, "Only Top", "lf");
        leg_store->AddEntry(HitRatio_Bottom_store, "Only Bottom", "lf");
        leg_store->Draw("same");


        // ========================================================
        // BAR PAD
        // ========================================================
        padBar_store->cd();

        TH1F* hFrame_store = (TH1F*)hFrame->Clone("hFrame_store");
        hFrame_store->Draw();

        double x1_store = x1;
        double x2_store = x2;
        double y0_store = 0.0;

        TLatex latex_store;
        latex_store.SetTextFont(132);
        latex_store.SetTextAlign(12);
        latex_store.SetTextSize(0.07);

        int smallLabelIndex_store = 0;

        auto drawFracBox_store = [&](double count, int color)
        {
            double frac = 100.0 * count / total;

            TBox* b = new TBox(x1_store, y0_store,
                            x2_store, y0_store + frac);
            b->SetFillColor(color);
            b->SetLineColor(kBlack);
            b->Draw("same");

            double yText = y0_store + 0.5 * frac;

            if (frac < 3.0) {
                yText = y0_store + frac + 3.0 * smallLabelIndex_store;
                smallLabelIndex_store++;
            }

            latex_store.SetTextColor(color);
            latex_store.DrawLatex(x2_store + 0.06,
                                yText,
                                Form("%.2f%%", frac));

            y0_store += frac;
        };

        // standard ROOT-safe colors
        drawFracBox_store(Full_count_bottom, kRed+1);
        drawFracBox_store(Full_count_top, kGreen+2);
        drawFracBox_store(Full_count_side, kBlue+1);
        drawFracBox_store(Full_count_top_and_bottom, kYellow+1);
        drawFracBox_store(Full_count_top_and_side, kCyan+1);
        drawFracBox_store(Full_count_bottom_and_side, kMagenta+1);

        gPad->RedrawAxis();

    // --------------------------------------------------------
    // Also store canvas inside plots_Reaction.root
    // --------------------------------------------------------

    TString rootFileName = Form("./%s_results/plots_%s.root",
                                ReactionName.c_str(),
                                ReactionName.c_str());

    TFile* plotFile = TFile::Open(rootFileName, "UPDATE");

    if (plotFile && !plotFile->IsZombie()) {

        // Create folder only if missing
        TDirectory* fissionDir =
            dynamic_cast<TDirectory*>(plotFile->Get("fission_plots"));

        if (!fissionDir) {
            fissionDir = plotFile->mkdir("fission_plots");
        }

        if (fissionDir) {
            fissionDir->cd();

            // overwrite old canvas if already there
            FissionEfficiencies->Write("FissionEfficiencies", TObject::kOverwrite);

            // Save image snapshot inside ROOT file
            TImage* img = TImage::Create();
            img->FromPad(cFission);
            img->Write("imgFission", TObject::kOverwrite);
            delete img;

            std::cout << "Saved FissionEfficiencies in "
                    << rootFileName
                    << "/fission_plots"
                    << std::endl;
        }

        plotFile->Close();
        delete plotFile;
    }
    else {
        std::cerr << "Could not open " << rootFileName
                << " in UPDATE mode." << std::endl;
    }

   
}