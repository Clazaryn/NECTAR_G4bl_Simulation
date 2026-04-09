#ifndef LUCAS_DRAWING_PALETTE_H
#define LUCAS_DRAWING_PALETTE_H

#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TAxis.h>
#include <TString.h>

namespace LucasStyle{

    static int DarkerBlueColor  = TColor::GetColor(53,39,140);
    static int DarkerBlueColor2 = TColor::GetColor(42,31,111);
    static int DarkerBlueColor3 = TColor::GetColor(31,23,83);
    static int BlueColor        = TColor::GetColor(41,95,166);
    static int LighterBlueColor = TColor::GetColor(56,128,223);

    static int RedColor         = TColor::GetColor(217,48,62);
    static int DarkerRedColor   = TColor::GetColor(129,28,37);
    static int BrownColor       = TColor::GetColor(72,16,21);

    static int GreenColor       = TColor::GetColor(10,166,122);
    static int DarkerGreenColor = TColor::GetColor(5,78,57);
    static int LightGreenColor  = TColor::GetColor(13,224,165);
    static int WaterGreenColor  = TColor::GetColor(0,224,138);

    static int OrangeColor      = TColor::GetColor(251,131,40);
    static int DarkerOrangeColor= TColor::GetColor(192,100,31);

    static int CyanColor        = TColor::GetColor(41,135,166);
    static int DarkerCyanColor  = TColor::GetColor(27,87,107);
    static int LighterCyanColor = TColor::GetColor(54,176,217);

    static int MagentaColor     = TColor::GetColor(218,48,164);
    static int DarkerMagentaColor = TColor::GetColor(129,28,97);

    static int YellowColor      = TColor::GetColor(237,196,59);
    static int PurpleColor      = TColor::GetColor(146,10,65);

    inline void SetGlobalStyle(){
        gStyle->SetPalette(57);
        gStyle->SetTextFont(132);
        gStyle->SetTitleFont(132,"xyz");
        gStyle->SetTitleFont(132,"a");
        gStyle->SetLabelFont(132,"xyz");
        gStyle->SetOptStat(0);
        gStyle->SetGridColor(kGray+2);
        gStyle->SetTitleFontSize(0.07);
        gROOT->ForceStyle();
    }

    inline void FormatHistSolid(TH1D* h, int col){
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        h->SetFillStyle(1001);
        h->SetFillColor(col);
        h->SetStats(0);
    }

    inline TGraphErrors* BuildDumbGraph(const char* graphTitle, const char* xTitle, const char* yTitle, double xmin, double xmax, double ymin, double ymax){
        TGraphErrors* g = new TGraphErrors();
        g->SetPoint(0, xmin, ymin);
        g->SetPoint(1, xmax, ymax);
        g->SetTitle(graphTitle);
        g->SetMarkerColor(0);
        g->SetMarkerStyle(8);

        g->GetXaxis()->SetTitle(xTitle);
        g->GetYaxis()->SetTitle(yTitle);

        g->GetXaxis()->CenterTitle(true);
        g->GetYaxis()->CenterTitle(true);

        g->GetXaxis()->SetRangeUser(xmin, xmax);
        g->GetYaxis()->SetRangeUser(ymin, ymax);

        g->GetXaxis()->SetTitleOffset(0.8);
        g->GetYaxis()->SetTitleOffset(0.8);

        g->GetXaxis()->SetTitleSize(0.055);
        g->GetYaxis()->SetTitleSize(0.055);

        g->GetXaxis()->SetLabelSize(0.045);
        g->GetYaxis()->SetLabelSize(0.045);

        return g;
    }

    inline TPad* BuildPad(const char* name, double x1, double y1, double x2, double y2, double leftMargin, double rightMargin, double bottomMargin, double topMargin){
        TPad* pad = new TPad(name, name, x1, y1, x2, y2);
        pad->SetLeftMargin(leftMargin);
        pad->SetRightMargin(rightMargin);
        pad->SetBottomMargin(bottomMargin);
        pad->SetTopMargin(topMargin);
        return pad;
    }

}

#endif