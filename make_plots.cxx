// --
// NECTAR plotting script
// Creates banana plots from det_analysis output files
// Last modification : G. Leckenby - 2025
//
// This file is compiled with g++ into a shared library

#include "UtilityScripts/make_plots.h"
#include "UtilityScripts/det_analysis.h"
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>

void make_plots() {
    // Load reaction information
    ReactionInfo reactionInfo;
    if (!reactionInfo.loadFromIni("reac_info.txt")) {
        std::cerr << "make_plots.cxx:" << __LINE__ << ": Error: Could not load reaction info. Exiting." << std::endl;
        return;
    }
    
    const char* reaction = reactionInfo.reaction.c_str();
    std::string det_setup = reactionInfo.det_setup;
    
    // Create one chain with ALL files (all channels, all excitation energies)
    TChain* fChain = new TChain("events");
    TString pattern = Form("./%s_results/Det_analysis/events_%s_*.root", reaction, reaction);
    Int_t nAdded = fChain->Add(pattern.Data());
    
    if (nAdded == 0) {
        std::cerr << "make_plots.cxx:" << __LINE__ << ": Error: No det_analysis files found for pattern: " << pattern << std::endl;
        delete fChain;
        return;
    }
    
    std::cout << "Added " << nAdded << " file(s) to chain. Note all ROOT files in Det_analysis directory will be added to the chain." << std::endl;
    
    // Create plot manager (no recType needed - it will handle channels internally)
    BananaPlotManager bananaPlots(fChain, det_setup, reaction, "");
    bananaPlots.initializePlots();
    bananaPlots.fillPlots();
    
    // Write plots (plot manager handles channel separation internally)
    bananaPlots.writePlots(reaction);
    
    delete fChain;
    
    std::cout << "\n=== All plotting complete ===" << std::endl;
}
