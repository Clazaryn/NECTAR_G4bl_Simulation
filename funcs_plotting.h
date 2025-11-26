#ifndef FUNCS_PLOTTING_H
#define FUNCS_PLOTTING_H

// ========= C++ libraries =========
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

// ========= ROOT libraries =========
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TCollection.h>
#include <TObject.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TBox.h>
#include <TLine.h>
#include <TMarker.h>

// Utility: Set X and Y axis ranges for first TH2 on a canvas
void SetTH2RangeOnCanvas(TCanvas* canvas, double xmin, double xmax, double ymin, double ymax);

// Create a custom color palette from two RGB colors
std::vector<int> CreateCustomPalette(const std::vector<double>& energies,
    double r1, double g1, double b1,  // Start color
    double r2, double g2, double b2   // End color
);

// Map excitation energies to color indices
std::map<double, int> MapEnergiesToColors(const std::vector<double>& energies,
    double r1, double g1, double b1,
    double r2, double g2, double b2
);

// Fit Gaussians and generate resolution graphs
TGraphErrors* FitResolutionGraph(TFile* file, const TString& prefix, TH2D* resolutionHist = nullptr);

// Fill accuracy histograms for reconstructed Eprot energy and theta
void FillAccuracyHistograms(TFile* file, const TString& prefix, TH2D* accuracyHist, const TString& suffix);

// Heavy residue coincidence plots in the XY plane
TH2D* histChannelPlots(const TString& reaction, const TString& prefix, const TString& suffix, TCanvas* canvas);

// Heavy residue coincidence survival fraction plots
void hrFractionPlot(TCanvas* canvas, const TString& prefix);

// escaping particles from the E1 and Eres detectors: veto plots for dE position map
void vetoPlots(TCanvas* canvas, const TString& prefix, const TString& reaction);

// Function to create veto-gated banana plots overlaid on ungated bananas
void vetoBananaPlots(TCanvas* canvasE1, TCanvas* canvasEres, const TString& prefix, const TString& reaction);

#endif // FUNCS_PLOTTING_H 