#ifndef MAKE_PLOTS_H
#define MAKE_PLOTS_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBox.h>
#include <TMultiGraph.h>

#include "reaction_info.h"
#include "det_analysis.h"

// ========= Base Plot Manager Class =========

class PlotManager {
protected:
    TChain* fChain;              // Chain of events to process
    std::string det_setup;       // "new" or "PoP"
    std::string reaction;
    std::string recType;
    
public:
    PlotManager(TChain* chain, const std::string& setup, 
                const std::string& reac, const std::string& rec);
    virtual ~PlotManager();
    
    // Initialize all histograms (must be implemented by derived classes)
    virtual void initializePlots() = 0;
    
    // Fill plots from event data (loops over chain, calls fillEvent for each)
    void fillPlots();
    
    // Fill a single event (must be implemented by derived classes)
    virtual void fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) = 0;
    
    // Write all plots to file (must be implemented by derived classes)
    // Takes reaction name - derived classes handle channel separation internally
    virtual void writePlots(const char* reaction) = 0;
};

// ========= Banana Plot Manager Class =========

class BananaPlotManager : public PlotManager {
private:
    // Histogram containers
    std::vector<TH2D*> hDEvE1_primary;      // dE vs E1 for primary telescope (New setup)
    std::vector<TH2D*> hDEvE1_auxillary;    // dE vs E1 for auxillary telescope (New setup)
    std::vector<TH2D*> hDEvE1_PoP;          // dE vs E1 for PoP setup
    
    std::vector<TH2D*> hDEvEres_primary;     // dE vs E_resid (punch-through), E_resid = E1+Eres (New)
    std::vector<TH2D*> hDEvEres_auxillary;   // dE vs E_resid for auxillary telescope (New setup)
    std::vector<TH2D*> hDEvEres_PoP;         // dE vs E_resid for PoP (E_resid = E1+E2+...)
    
    // Theta-binned histograms (New setup only)
    std::vector<std::vector<TH2D*>> hDEvE1_primary_theta;     // [theta_bin]
    std::vector<std::vector<TH2D*>> hDEvE1_auxillary_theta;   // [theta_bin]
    std::vector<std::vector<TH2D*>> hDEvEres_primary_theta;   // [theta_bin]
    std::vector<std::vector<TH2D*>> hDEvEres_auxillary_theta; // [theta_bin]
    
    // Theta binning: primary 35-45,...,75-85 deg; auxillary 5-15, 15-25 deg (avoid empty bins)
    static const Int_t nThetaBins_primary = 5;    // 35-45, 45-55, 55-65, 65-75, 75-85
    static const Int_t nThetaBins_auxillary = 2;  // 5-15, 15-25
    Double_t thetaBinEdges_primary[nThetaBins_primary + 1];
    Double_t thetaBinEdges_auxillary[nThetaBins_auxillary + 1];
    
    // Histogram parameters
    Int_t nBins_dE, nBins_E1, nBins_Eres, nBins_Eresid;
    Double_t dE_min, dE_max, E1_min, E1_max, Eres_min, Eres_max, Eresid_max;
    
    Int_t getThetaBinPrimary(Double_t theta) const;
    Int_t getThetaBinAuxillary(Double_t theta) const;
    
public:
    BananaPlotManager(TChain* chain, const std::string& setup, 
                     const std::string& reac, const std::string& rec);
    virtual ~BananaPlotManager();
    
    // Initialize all histograms
    void initializePlots() override;
    
    // Fill a single event
    void fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) override;
    
    // Write all plots to file
    void writePlots(const char* reaction) override;
};

// ========= Accuracy Plot Manager Class =========
// Compares reconstructed vs true values: Eexc, Eejc (ejectile total energy), theta.
// New setup: separate primary and auxillary (routed by detector_id); subfolders primary/ and auxillary/.
// PoP setup: single set in PoP/ subfolder.
// Per detector: 3x 1D accuracy + 9x 2D accuracy vs E*, Eejc, theta.

class AccuracyPlotManager : public PlotManager {
private:
    // Index convention: 1D [0]=dEexc, [1]=dEejc, [2]=dTheta
    // 2D [0-2]=vs Eexc, [3-5]=vs Eejc, [6-8]=vs Theta (each triple: dEexc, dEejc, dTheta)
    std::vector<TH1D*> h_1d_primary;    // New setup: primary telescope (theta 35-85 deg)
    std::vector<TH2D*> h_2d_primary;
    std::vector<TH1D*> h_1d_auxillary;  // New setup: auxillary telescope (theta 5-25 deg)
    std::vector<TH2D*> h_2d_auxillary;
    std::vector<TH1D*> h_1d_PoP;         // PoP setup (theta 51-69 deg)
    std::vector<TH2D*> h_2d_PoP;

public:
    AccuracyPlotManager(TChain* chain, const std::string& setup,
                       const std::string& reac, const std::string& rec);
    virtual ~AccuracyPlotManager();

    void initializePlots() override;
    void fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) override;
    void writePlots(const char* reaction) override;
};

// ========= Excitation Energy Resolution Plot Manager Class =========
// 2D histogram of (recon_Eexc - true_Eexc) vs true_Eexc; fit Gaussian per E* slice to get sigma(E*).
// Output: resolution map (2D) and TGraphErrors sigma vs E* in exc_resolution/ folder.

class ExcResolutionPlotManager : public PlotManager {
private:
    TH2D* h_dEexc_vs_Eexc_primary;    // New: primary telescope
    TH2D* h_dEexc_vs_Eexc_auxillary;  // New: auxillary telescope
    TH2D* h_dEexc_vs_Eexc_PoP;        // PoP: single telescope
    TGraphErrors* g_resolution_primary;
    TGraphErrors* g_resolution_auxillary;
    TGraphErrors* g_resolution_PoP;
    TCanvas* c_resolution;

public:
    ExcResolutionPlotManager(TChain* chain, const std::string& setup,
                            const std::string& reac, const std::string& rec);
    virtual ~ExcResolutionPlotManager();

    void initializePlots() override;
    void fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) override;
    void writePlots(const char* reaction) override;
};

// ========= Transmission Plot Manager Class =========
// Coincidence 2D scatter (MagSept, HRplane, QuadWall) and transmission fraction vs E*.
// New: primary/auxillary folders. PoP: single folder. Channels: gamma(red), 1n(blue), 2n(green), 3n(purple), 4n(cyan).

class TransmissionPlotManager : public PlotManager {
private:
    // Coincidence: [virt][channel][Eexc_bin] TH2D. virt: 0=MagSept, 1=HRplane, 2=QuadWall
    static const Int_t nChannels = 5;
    static const Int_t nEexcBins = 5;
    std::vector<std::vector<std::vector<TH2D*>>> h_coinc_primary;
    std::vector<std::vector<std::vector<TH2D*>>> h_coinc_auxillary;
    std::vector<std::vector<std::vector<TH2D*>>> h_coinc_PoP;
    // Transmission: TH1D binned in Eexc per channel for MagSept and HRplane
    std::vector<TH1D*> h_MagSept_primary, h_HRplane_primary;
    std::vector<TH1D*> h_MagSept_auxillary, h_HRplane_auxillary;
    std::vector<TH1D*> h_MagSept_PoP, h_HRplane_PoP;
    TCanvas* c_MagSept_primary, * c_HRplane_primary, * c_QuadWall_primary;
    TCanvas* c_MagSept_auxillary, * c_HRplane_auxillary, * c_QuadWall_auxillary;
    TCanvas* c_MagSept_PoP, * c_HRplane_PoP, * c_QuadWall_PoP;
    TCanvas* c_transmission_primary, * c_transmission_auxillary, * c_transmission_PoP;

public:
    TransmissionPlotManager(TChain* chain, const std::string& setup,
                           const std::string& reac, const std::string& rec);
    virtual ~TransmissionPlotManager();

    void initializePlots() override;
    void fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) override;
    void writePlots(const char* reaction) override;
};

// ========= Utility Functions =========

// Detect telescope setup from reac_info.txt
std::string detectTelescopeSetup(const std::string& reac_info_file = "reac_info.txt");

// ========= Main plotting function =========
void make_plots();

#endif // MAKE_PLOTS_H
