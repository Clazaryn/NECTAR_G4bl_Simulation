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
#include <TH2D.h>
#include <TDirectory.h>

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
    
    std::vector<TH2D*> hDEvEres_primary;     // dE vs Eres for primary telescope (New setup)
    std::vector<TH2D*> hDEvEres_auxillary;  // dE vs Eres for auxillary telescope (New setup)
    std::vector<TH2D*> hDEvEres_PoP;        // dE vs Eres for PoP setup
    
    // Theta-binned histograms (New setup only)
    std::vector<std::vector<TH2D*>> hDEvE1_primary_theta;     // [theta_bin]
    std::vector<std::vector<TH2D*>> hDEvE1_auxillary_theta;   // [theta_bin]
    std::vector<std::vector<TH2D*>> hDEvEres_primary_theta;     // [theta_bin]
    std::vector<std::vector<TH2D*>> hDEvEres_auxillary_theta; // [theta_bin]
    
    // Theta binning parameters
    static const Int_t nThetaBins = 18;  // 0-180 degrees in 10 degree bins
    Double_t thetaBinEdges[nThetaBins + 1];
    
    // Histogram parameters
    Int_t nBins_dE, nBins_E1, nBins_Eres;
    Double_t dE_min, dE_max, E1_min, E1_max, Eres_min, Eres_max;
    
    // Helper function to get theta bin index
    Int_t getThetaBin(Double_t theta) const;
    
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

// ========= Utility Functions =========

// Detect telescope setup from reac_info.txt
std::string detectTelescopeSetup(const std::string& reac_info_file = "reac_info.txt");

// ========= Main plotting function =========
void make_plots();

#endif // MAKE_PLOTS_H
