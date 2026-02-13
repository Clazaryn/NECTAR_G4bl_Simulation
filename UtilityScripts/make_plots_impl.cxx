// Implementation file for make_plots.h (base class)
#include "make_plots.h"
#include <iostream>

// ========= PlotManager Base Class Implementation =========

PlotManager::PlotManager(TChain* chain, const std::string& setup, 
                         const std::string& reac, const std::string& rec)
    : fChain(chain), det_setup(setup), reaction(reac), recType(rec) {
}

PlotManager::~PlotManager() {
    // Don't delete fChain here - it's managed externally
}

void PlotManager::fillPlots() {
    if (!fChain) {
        std::cerr << "make_plots_impl.cxx:" << __LINE__ << ": Error: fChain is null" << std::endl;
        return;
    }
    
    // Set up branches
    LightEjectile* ejectile = new LightEjectile();
    HeavyResidue* residue = new HeavyResidue();
    Short_t decay_channel = -1;
    
    fChain->SetBranchAddress("ejectile", &ejectile);
    fChain->SetBranchAddress("residue", &residue);
    fChain->SetBranchAddress("decay_channel", &decay_channel);
    
    // Loop over all events
    Long64_t nEntries = fChain->GetEntries();
    std::cout << "Processing " << nEntries << " events from chain..." << std::endl;
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        fChain->GetEntry(i);
        fillEvent(ejectile, residue, decay_channel);
    }
    
    delete ejectile;
    delete residue;
    std::cout << "Finished filling plots" << std::endl;
}

// ========= Utility Functions =========

std::string detectTelescopeSetup(const std::string& reac_info_file) {
    ReactionInfo reactionInfo;
    if (!reactionInfo.loadFromIni(reac_info_file)) {
        std::cerr << "make_plots_impl.cxx:" << __LINE__ << ": Error: Could not load reaction info from " << reac_info_file << std::endl;
        return "";
    }
    return reactionInfo.det_setup;
}
