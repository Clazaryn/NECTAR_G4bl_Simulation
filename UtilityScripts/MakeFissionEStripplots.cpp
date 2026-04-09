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


#include "Lucas_drawing_palette.h"

using namespace std;
//namespace fs = std::filesystem;

double* Excitationenergies;
double EexcBinSize;
int Bin_number;



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

void FillFragmentPlots(TH2D* E_depositedVsTheta[4][3], TH2D* E_depositedVsStripVert[4][3], TH2D* E_depositedVsStripHor[4][3], TH2D* Heatmap[4][3], int detIdx, int fragIdx, const FissionFragment& frag){
    
    //cout << "___test__" << frag.true_theta << endl;

    E_depositedVsTheta[detIdx][fragIdx]->Fill(frag.true_theta, frag.E_deposited);
    E_depositedVsStripVert[detIdx][fragIdx]->Fill(frag.vert_strip, frag.E_deposited);
    E_depositedVsStripHor[detIdx][fragIdx]->Fill(frag.hor_strip, frag.E_deposited);
    Heatmap[detIdx][fragIdx]->Fill(frag.vert_strip, frag.hor_strip);
}

//____________________________________________________________

void ProcessFragmentForDetector(bool hit, int detIdx, int fragIdx, const FissionFragment& frag, TH2D* E_depositedVsTheta[4][3], TH2D* E_depositedVsStripVert[4][3], TH2D* E_depositedVsStripHor[4][3], TH2D* Heatmap[4][3]){
    if(!hit) return;
    // specific detector + specific fragment 
    FillFragmentPlots(E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap, detIdx, fragIdx, frag);
    // specific detector + all fragments
    FillFragmentPlots(E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap, detIdx, 2, frag);
    // all detectors + specific fragment
    FillFragmentPlots(E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap, 3, fragIdx, frag);
    // all detectors + all fragments
    FillFragmentPlots(E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap, 3, 2, frag);
}

//___Helper__for_this_specific_drawing______________________________

const char* GetFragmentTypeName(int type){
    if(type==0) return "Heavy fragments";
    if(type==1) return "Light fragments";
    return "All fragments";
}

const char* GetDetectorName(int det){
    if(det==0) return "Bottom detector";
    if(det==1) return "Side detector";
    if(det==2) return "Top detector";
    return "All detectors";
}

//___Helper__for_Graphs_and_Histo______________________________

void DrawEdepThetaPad(TPad* pad, TH2D* hist, const char* detName, double xmin, double xmax, double ymin, double ymax){

    TGraphErrors* DumbGraph = LucasStyle::BuildDumbGraph("", "#Theta (Deg.)", "E_{dep} (MeV)", xmin, xmax, ymin, ymax);
    DumbGraph->Draw("AP");

    hist->Draw("COL SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(132);
    latex.SetTextSize(0.060);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.45, 0.98, detName);

    pad->SetGrid();
    gPad->RedrawAxis("G");
    gPad->RedrawAxis();
}

void DrawEdepHorStripPad(TPad* pad, TH2D* hist, const char* detName, double xmin, double xmax, double ymin, double ymax){

    TGraphErrors* DumbGraph = LucasStyle::BuildDumbGraph("", "Horizontal Strip Number", "E_{dep} (MeV)", xmin, xmax, ymin, ymax);
    DumbGraph->Draw("AP");

    hist->RebinY(2);
    hist->Draw("COL SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(132);
    latex.SetTextSize(0.060);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.45, 0.98, detName);

    pad->SetGrid();
    gPad->RedrawAxis("G");
    gPad->RedrawAxis();
}

void DrawEdepVertStripPad(TPad* pad, TH2D* hist, const char* detName, double xmin, double xmax, double ymin, double ymax){

    TGraphErrors* DumbGraph = LucasStyle::BuildDumbGraph("", "Vertical Strip Number", "E_{dep} (MeV)", xmin, xmax, ymin, ymax);
    DumbGraph->Draw("AP");

    hist->RebinY(2);
    hist->Draw("COL SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(132);
    latex.SetTextSize(0.060);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.45, 0.98, detName);

    pad->SetGrid();
    gPad->RedrawAxis("G");
    gPad->RedrawAxis();
}

void DrawHeatMapPad(TPad* pad, TH2D* hist, const char* detName, double xmin, double xmax, double ymin, double ymax){

    TGraphErrors* DumbGraph = LucasStyle::BuildDumbGraph("", "Vertical Strip Number", "Horizontal Strip Number", xmin, xmax, ymin, ymax);
    DumbGraph->Draw("AP");

    hist->Draw("COL SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(132);
    latex.SetTextSize(0.060);
    latex.SetTextAlign(13);
    latex.DrawLatex(0.45, 0.98, detName);

    pad->SetGrid();
    gPad->RedrawAxis("G");
    gPad->RedrawAxis();
}

//___Helper__for_canvas______________________________


TCanvas* BuildEdepThetaCanvas(const char* canvasName, const char* canvasTitle, TH2D* E_depositedVsTheta[4][3], int type){
    double xmin = 0.0;
    double xmax = 19.0;
    double ymin = 0.0;
    double ymax = 3000.0;

    double titleY1 = 0.92;

    double leftMargin = 0.12;
    double rightMargin = 0.06;
    double bottomMargin = 0.12;
    double topMargin = 0.10;

    TCanvas* c = new TCanvas(canvasName, canvasTitle, 1920, 1080);

    TPad* padTitle  = LucasStyle::BuildPad(Form("%s_padTitle",canvasName),  0.00, titleY1, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00);
    TPad* padTop    = LucasStyle::BuildPad(Form("%s_padTop",canvasName),    0.00, 0.46, 0.50, titleY1, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padSide   = LucasStyle::BuildPad(Form("%s_padSide",canvasName),   0.50, 0.46, 1.00, titleY1, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padBottom = LucasStyle::BuildPad(Form("%s_padBottom",canvasName), 0.00, 0.00, 0.50, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padAll    = LucasStyle::BuildPad(Form("%s_padAll",canvasName),    0.50, 0.00, 1.00, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);

    padTitle->Draw();
    padTop->Draw();
    padSide->Draw();
    padBottom->Draw();
    padAll->Draw();

    padTitle->cd();

    TLatex bigTitle;
    bigTitle.SetNDC();
    bigTitle.SetTextFont(132);
    bigTitle.SetTextSize(0.45);
    bigTitle.SetTextAlign(22);
    bigTitle.DrawLatex(0.50, 0.45, Form("E_{dep} as a function of #Theta emission angle - %s", GetFragmentTypeName(type)));

    padTop->cd();
    DrawEdepThetaPad(padTop,    E_depositedVsTheta[2][type], "Top detector",    xmin, xmax, ymin, ymax);
    padSide->cd();
    DrawEdepThetaPad(padSide,   E_depositedVsTheta[1][type], "Side detector",   xmin, xmax, ymin, ymax);
    padBottom->cd();
    DrawEdepThetaPad(padBottom, E_depositedVsTheta[0][type], "Bottom detector", xmin, xmax, ymin, ymax);
    padAll->cd();
    DrawEdepThetaPad(padAll,    E_depositedVsTheta[3][type], "All detectors",   xmin, xmax, ymin, ymax);

    c->Update();
    return c;
}

TCanvas* BuildEdepHorStripCanvas(const char* canvasName, const char* canvasTitle, TH2D* E_depositedVsTheta[4][3], int type){
    double xmin = 0.0;
    double xmax = 17.0;
    double xmax_side = 41.0;
    double ymin = 0.0;
    double ymax = 3000.0;

    double titleY1 = 0.92;

    double leftMargin = 0.12;
    double rightMargin = 0.06;
    double bottomMargin = 0.12;
    double topMargin = 0.10;

    TCanvas* c = new TCanvas(canvasName, canvasTitle, 1920, 1080);

    TPad* padTitle  = LucasStyle::BuildPad(Form("%s_padTitle",canvasName),  0.00, titleY1, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00);
    TPad* padTop    = LucasStyle::BuildPad(Form("%s_padTop",canvasName),    0.00, 0.46, 0.50, titleY1, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padSide   = LucasStyle::BuildPad(Form("%s_padSide",canvasName),   0.50, 0.26, 1.00, 0.72, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padBottom = LucasStyle::BuildPad(Form("%s_padBottom",canvasName), 0.00, 0.00, 0.50, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);
    //TPad* padAll    = LucasStyle::BuildPad(Form("%s_padAll",canvasName),    0.50, 0.00, 1.00, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);

    padTitle->Draw();
    padTop->Draw();
    padSide->Draw();
    padBottom->Draw();
    //padAll->Draw();

    padTitle->cd();

    TLatex bigTitle;
    bigTitle.SetNDC();
    bigTitle.SetTextFont(132);
    bigTitle.SetTextSize(0.45);
    bigTitle.SetTextAlign(22);
    bigTitle.DrawLatex(0.50, 0.45, Form("E_{dep} as a function of Horizontal strip number - %s", GetFragmentTypeName(type)));

    padTop->cd();
    DrawEdepHorStripPad(padTop,    E_depositedVsTheta[2][type], "Top detector",    xmin, xmax, ymin, ymax);
    padSide->cd();
    DrawEdepHorStripPad(padSide,   E_depositedVsTheta[1][type], "Side detector",   xmin, xmax_side, ymin, ymax);
    padBottom->cd();
    DrawEdepHorStripPad(padBottom, E_depositedVsTheta[0][type], "Bottom detector", xmin, xmax, ymin, ymax);
    //padAll->cd();
    //DrawEdepThetaPad(padAll,    E_depositedVsTheta[3][type], "All detectors",   xmin, xmax, ymin, ymax);

    c->Update();
    return c;
}

TCanvas* BuildEdepVertStripCanvas(const char* canvasName, const char* canvasTitle, TH2D* E_depositedVsTheta[4][3], int type){
    double xmin = 0.0;
    double xmax = 17.0;
    double xmax_side = 81.0;
    double ymin = 0.0;
    double ymax = 3000.0;

    double titleY1 = 0.92;

    double leftMargin = 0.12;
    double rightMargin = 0.06;
    double bottomMargin = 0.12;
    double topMargin = 0.10;

    TCanvas* c = new TCanvas(canvasName, canvasTitle, 1920, 1080);

    TPad* padTitle  = LucasStyle::BuildPad(Form("%s_padTitle",canvasName),  0.00, titleY1, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00);
    TPad* padTop    = LucasStyle::BuildPad(Form("%s_padTop",canvasName),    0.00, 0.46, 0.50, titleY1, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padSide   = LucasStyle::BuildPad(Form("%s_padSide",canvasName),   0.50, 0.26, 1.00, 0.72, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padBottom = LucasStyle::BuildPad(Form("%s_padBottom",canvasName), 0.00, 0.00, 0.50, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);
    //TPad* padAll    = LucasStyle::BuildPad(Form("%s_padAll",canvasName),    0.50, 0.00, 1.00, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);

    padTitle->Draw();
    padTop->Draw();
    padSide->Draw();
    padBottom->Draw();
    //padAll->Draw();

    padTitle->cd();

    TLatex bigTitle;
    bigTitle.SetNDC();
    bigTitle.SetTextFont(132);
    bigTitle.SetTextSize(0.45);
    bigTitle.SetTextAlign(22);
    bigTitle.DrawLatex(0.50, 0.45, Form("E_{dep} as a function of Vertical strip number - %s", GetFragmentTypeName(type)));

    padTop->cd();
    DrawEdepVertStripPad(padTop,    E_depositedVsTheta[2][type], "Top detector",    xmin, xmax, ymin, ymax);
    padSide->cd();
    DrawEdepVertStripPad(padSide,   E_depositedVsTheta[1][type], "Side detector",   xmin, xmax_side, ymin, ymax);
    padBottom->cd();
    DrawEdepVertStripPad(padBottom, E_depositedVsTheta[0][type], "Bottom detector", xmin, xmax, ymin, ymax);
    //padAll->cd();
    //DrawEdepThetaPad(padAll,    E_depositedVsTheta[3][type], "All detectors",   xmin, xmax, ymin, ymax);

    c->Update();
    return c;
}

TCanvas* BuildHeatMapCanvas(const char* canvasName, const char* canvasTitle, TH2D* E_depositedVsTheta[4][3], int type){
    double xmin = 0.0;
    double xmax = 17.0;
    double xmax_side = 81.0;
    double ymin = 0.0;
    double ymax = 17.0;
    double ymax_side = 45.0;

    double titleY1 = 0.92;

    double leftMargin = 0.12;
    double rightMargin = 0.06;
    double bottomMargin = 0.12;
    double topMargin = 0.10;

    TCanvas* c = new TCanvas(canvasName, canvasTitle, 1920, 1080);

    TPad* padTitle  = LucasStyle::BuildPad(Form("%s_padTitle",canvasName),  0.00, titleY1, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00);
    TPad* padTop    = LucasStyle::BuildPad(Form("%s_padTop",canvasName),    0.00, 0.46, 0.50, titleY1, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padSide   = LucasStyle::BuildPad(Form("%s_padSide",canvasName),   0.50, 0.26, 1.00, 0.72, leftMargin, rightMargin, bottomMargin, topMargin);
    TPad* padBottom = LucasStyle::BuildPad(Form("%s_padBottom",canvasName), 0.00, 0.00, 0.50, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);
    //TPad* padAll    = LucasStyle::BuildPad(Form("%s_padAll",canvasName),    0.50, 0.00, 1.00, 0.46, leftMargin, rightMargin, bottomMargin, topMargin);

    padTitle->Draw();
    padTop->Draw();
    padSide->Draw();
    padBottom->Draw();
    //padAll->Draw();

    padTitle->cd();

    TLatex bigTitle;
    bigTitle.SetNDC();
    bigTitle.SetTextFont(132);
    bigTitle.SetTextSize(0.45);
    bigTitle.SetTextAlign(22);
    bigTitle.DrawLatex(0.50, 0.45, Form("Heat Maps - %s", GetFragmentTypeName(type)));

    padTop->cd();
    DrawHeatMapPad(padTop,    E_depositedVsTheta[2][type], "Top detector",    xmin, xmax, ymin, ymax);
    padSide->cd();
    DrawHeatMapPad(padSide,   E_depositedVsTheta[1][type], "Side detector",   xmin, xmax_side, ymin, ymax_side);
    padBottom->cd();
    DrawHeatMapPad(padBottom, E_depositedVsTheta[0][type], "Bottom detector", xmin, xmax, ymin, ymax);
    //padAll->cd();
    //DrawEdepThetaPad(padAll,    E_depositedVsTheta[3][type], "All detectors",   xmin, xmax, ymin, ymax);

    c->Update();
    return c;
}

//__Helper to save canvas_____

void SaveCanvasInRootAndPDF(TCanvas* c, const string& ReactionName, const char* canvasWriteName){
    TString rootFileName = Form("./%s_results/plots_%s.root", ReactionName.c_str(), ReactionName.c_str());
    TString pdfDirName = Form("./%s_results/Fission_pdfs", ReactionName.c_str());
    TString pdfFileName = Form("./%s_results/Fission_pdfs/%s.pdf", ReactionName.c_str(), canvasWriteName);

    if(gSystem->AccessPathName(pdfDirName)){
        gSystem->mkdir(pdfDirName, kTRUE);
    }

    TFile* plotFile = TFile::Open(rootFileName, "UPDATE");

    if(plotFile && !plotFile->IsZombie()){
        TDirectory* fissionDir = dynamic_cast<TDirectory*>(plotFile->Get("fission_plots"));

        if(!fissionDir){
            fissionDir = plotFile->mkdir("fission_plots");
        }

        if(fissionDir){
            fissionDir->cd();
            c->Write(canvasWriteName, TObject::kOverwrite);

            std::cout << "Saved canvas " << canvasWriteName << " in " << rootFileName << "/fission_plots" << std::endl;
        }

        plotFile->Close();
        delete plotFile;
    }
    else{
        std::cerr << "Could not open " << rootFileName << " in UPDATE mode." << std::endl;
    }

    c->SaveAs(pdfFileName);
    std::cout << "Saved pdf " << pdfFileName << std::endl;
}

//____________________________________________________________
void MakeFissionEStripplots(const char* reaction_name = "238U_dp"){

    //string ReactionName = "238U_dp";
    string ReactionName = reaction_name;
    EexcBinSize = 1.0;

    // -------------------------------------------------------------------------
    // Let's find out what output file we have and loop through them 
    // -------------------------------------------------------------------------

    int energyNumber;
    Excitationenergies = nullptr;

    bool ok = FindFissionEnergies("./"+ReactionName+"_results/Det_analysis",
                                energyNumber,
                                Excitationenergies);

    energyNumber--;

    double HalfEnergyStep = (Excitationenergies[1]-Excitationenergies[0])/2.0;

    // --------------------------------------------------------
    // Set style
    // --------------------------------------------------------

    LucasStyle::SetGlobalStyle();

    // ---------------------------------------------------------------------------------------------------------------------
    // Let's loop on the fission ".root" files and gets the number of events (and proportion) for each of the detection case  
    // ---------------------------------------------------------------------------------------------------------------------
    

    LightEjectile* ejectile = new LightEjectile();
    FissionEvent* fission = new FissionEvent();

    //______For the ejectile____________

    Int_t Z, A;
    Int_t detector_id;  // 0=primary, 1=auxillary for New detectors
    Int_t vert_strip, hor_strip;
    Double_t true_Eexc, true_Eejc, true_theta;
    Double_t recon_Eexc, recon_Eejc, recon_theta;
    Float_t meas_dE, meas_E1, meas_Eres;  // Measured Eres energy (crystal in New setup)
    Float_t meas_E2, meas_E3, meas_E4, meas_E5, meas_E6;  // Measured E2 through E6 energy (thick Si in PoP setup)
    
    //________For the fission fragments_____________
    
    Double_t  detec_x, detec_y; 
    Double_t E_emission;
    Double_t E_deposited;

    Bool_t light_hit_Sidedetec;
    Bool_t light_hit_Topdetec;
    Bool_t light_hit_Bottomdetec;

    Bool_t heavy_hit_Sidedetec;
    Bool_t heavy_hit_Topdetec;
    Bool_t heavy_hit_Bottomdetec;

    Int_t light_Z, light_A;
    Int_t heavy_Z, heavy_A; 


    //cout << "__from__" << Excitationenergies[0] << "__to__" << Excitationenergies[energyNumber] << "___binSize__" << EexcBinSize << "___bin_number__" 
    //     << Bin_number << "__check__" <<  Excitationenergies[energyNumber] << "__==__" << Excitationenergies[0]+(EexcBinSize*Bin_number) << endl;

    // First coordinate : 
    // 0 = Bottom, 1 = Side, 2 = Top, 3 = All => Not really relevant for strip btw but who knows, for future symmetric geometry ?
    // Second coordinate : 
    // 0 = Heavy, 1 = Light, 2 = All 
    TH2D* E_depositedVsTheta[4][3];
    TH2D* E_depositedVsStripVert[4][3];
    TH2D* E_depositedVsStripHor[4][3];
    TH2D* Heatmap[4][3];


    for(int type=0; type<3;type++){
        for(int det=0; det<4;det++){

            TString nameEtheta = Form("EVsTheta_%02d_%02d",det,type);
            TString nameEStripVert = Form("EVsStripVert_%02d_%02d",det,type);
            TString nameEStripHor = Form("EVsStripHor_%02d_%02d",det,type);
            TString nameHeatmap = Form("Heatmap_%02d_%02d",det,type);

            E_depositedVsTheta[det][type] = new TH2D(nameEtheta, nameEtheta, 200, 0, 20, 150, 0, 3000);
            E_depositedVsStripVert[det][type] = new TH2D(nameEStripVert, nameEStripVert, 80, 0, 80, 150, 0, 3000);
            E_depositedVsStripHor[det][type] = new TH2D(nameEStripHor, nameEStripHor, 80, 0, 80, 150, 0, 3000);
            Heatmap[det][type] = new TH2D(nameHeatmap, nameHeatmap, 80, 0, 80, 80, 0, 80);
        }
    }
    

    
    // Let's find and open the tree associated with the reaction ! 
    TChain * tree = new TChain("events");

    for (int Eex_bin = 0; Eex_bin < energyNumber; Eex_bin++){
        if(Excitationenergies[Eex_bin]>9.999){
            tree->Add(("./"+ReactionName+"_results/Det_analysis/events_"+ReactionName+Form("_HRf_excEn%.1fMeV.root",Excitationenergies[Eex_bin])).c_str());
        }
        if(Excitationenergies[Eex_bin]<10){
            tree->Add(("./"+ReactionName+"_results/Det_analysis/events_"+ReactionName+Form("_HRf_excEn0%.1fMeV.root",Excitationenergies[Eex_bin])).c_str());
        }
    }

    tree->SetBranchAddress("ejectile", &ejectile);
    tree->SetBranchAddress("fission", &fission);

    Long64_t entries = tree->GetEntries();

    
    //Let's go through the tree now !!  
    for (Long64_t i = 0; i < entries; i++) {

        if (i % 5000 == 0) {
            double prog = 100.0 * (double)i / (double)entries;
            cout << "\rProgression " << prog << " %" << flush;
        }

        tree->GetEntry(i);  

        //if(ejectile->true_Eexc>0){
        //    cout << "__entry_" << i << "__ejectile_E*__" << ejectile->true_Eexc 
        //         << "__light__hit(Top,Bottom,Side)__(" << fission->light.hit_Topdetec << ", " << fission->light.hit_Bottomdetec << ", " << fission->light.hit_Sidedetec << ") " 
        //         << "__heavy__hit(Top,Bottom,Side)__(" << fission->heavy.hit_Topdetec << ", " << fission->heavy.hit_Bottomdetec << ", " << fission->heavy.hit_Sidedetec << ") " 
        //         << endl;
        //}


        //cout << "____test__" << fission->heavy.hit_Bottomdetec << endl;

        // Bottom detector = index 0
        ProcessFragmentForDetector(fission->heavy.hit_Bottomdetec, 0, 0, fission->heavy, E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap);
        ProcessFragmentForDetector(fission->light.hit_Bottomdetec, 0, 1, fission->light, E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap);

        // Side detector = index 1
        ProcessFragmentForDetector(fission->heavy.hit_Sidedetec, 1, 0, fission->heavy, E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap);
        ProcessFragmentForDetector(fission->light.hit_Sidedetec, 1, 1, fission->light, E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap);

        // Top detector = index 2
        ProcessFragmentForDetector(fission->heavy.hit_Topdetec, 2, 0, fission->heavy, E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap);
        ProcessFragmentForDetector(fission->light.hit_Topdetec, 2, 1, fission->light, E_depositedVsTheta, E_depositedVsStripVert, E_depositedVsStripHor, Heatmap);

    }

        
    delete tree;
    tree = nullptr;   // good practice

    //E_depositedVsTheta[0][0]->Draw("colz");
    TCanvas* cFissionEHorStrip_Heavy = BuildEdepHorStripCanvas("cFissionEHorStrip_Heavy", "Fission Edep vs HorStrip Heavy", E_depositedVsStripHor, 0);
    TCanvas* cFissionEHorStrip_Light = BuildEdepHorStripCanvas("cFissionEHorStrip_Light", "Fission Edep vs HorStrip Light", E_depositedVsStripHor, 1);
    TCanvas* cFissionEHorStrip_All = BuildEdepHorStripCanvas("cFissionEHorStrip_All", "Fission Edep vs HorStrip All", E_depositedVsStripHor, 2);

    TCanvas* cFissionEVertStrip_Heavy = BuildEdepVertStripCanvas("cFissionEVertStrip_Heavy", "Fission Edep vs VertStrip Heavy", E_depositedVsStripVert, 0);
    TCanvas* cFissionEVertStrip_Light = BuildEdepVertStripCanvas("cFissionEVertStrip_Light", "Fission Edep vs VertStrip Light", E_depositedVsStripVert, 1);
    TCanvas* cFissionEVertStrip_All = BuildEdepVertStripCanvas("cFissionEVertStrip_All", "Fission Edep vs VertStrip All", E_depositedVsStripVert, 2);

    TCanvas* cHeatmap_Heavy = BuildHeatMapCanvas("cFissionHeatMap_Heavy", "Fission Horstrip Vs VertStrip Heavy", Heatmap, 0);
    TCanvas* cHeatmap_Light = BuildHeatMapCanvas("cFissionHeatMap_Light", "Fission Horstrip Vs VertStrip Light", Heatmap, 1);
    TCanvas* cHeatmap_All = BuildHeatMapCanvas("cFissionHeatMap_All", "Fission Horstrip Vs VertStrip All", Heatmap, 2);

    TCanvas* cFissionETheta_Heavy = BuildEdepThetaCanvas("cFissionETheta_Heavy", "Fission Edep vs Theta Heavy", E_depositedVsTheta, 0);
    TCanvas* cFissionETheta_Light = BuildEdepThetaCanvas("cFissionETheta_Light", "Fission Edep vs Theta Light", E_depositedVsTheta, 1);
    TCanvas* cFissionETheta_All = BuildEdepThetaCanvas("cFissionETheta_All", "Fission Edep vs Theta All", E_depositedVsTheta, 2);

    SaveCanvasInRootAndPDF(cFissionETheta_Heavy, ReactionName, "Fission_ETheta_Heavy");
    SaveCanvasInRootAndPDF(cFissionETheta_Light, ReactionName, "Fission_ETheta_Light");
    SaveCanvasInRootAndPDF(cFissionETheta_All, ReactionName, "Fission_ETheta_All");

    SaveCanvasInRootAndPDF(cFissionEHorStrip_Heavy, ReactionName, "Fission_EHorStrip_Heavy");
    SaveCanvasInRootAndPDF(cFissionEHorStrip_Light, ReactionName, "Fission_EHorStrip_Light");
    SaveCanvasInRootAndPDF(cFissionEHorStrip_All, ReactionName, "Fission_EHorStrip_All");

    SaveCanvasInRootAndPDF(cFissionEVertStrip_Heavy, ReactionName, "Fission_EVertStrip_Heavy");
    SaveCanvasInRootAndPDF(cFissionEVertStrip_Light, ReactionName, "Fission_EVertStrip_Light");
    SaveCanvasInRootAndPDF(cFissionEVertStrip_All, ReactionName, "Fission_EVertStrip_All");

    SaveCanvasInRootAndPDF(cHeatmap_Heavy, ReactionName, "Fission_HeatMap_Heavy");
    SaveCanvasInRootAndPDF(cHeatmap_Light, ReactionName, "Fission_HeatMap_Light");
    SaveCanvasInRootAndPDF(cHeatmap_All, ReactionName, "Fission_HeatMap_All");
    
    


   
}