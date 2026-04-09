// Implementation file for make_plots.h (base class)
#include "make_plots.h"
#include "ini_parser.h"
#include <iostream>
#include <cmath>

// ========= PlotManager Base Class Implementation =========

PlotManager::PlotManager(TChain* chain, const std::string& setup, 
                         const std::string& reac, const std::string& rec)
    : fChain(chain), det_setup(setup), reaction(reac), recType(rec) {
}

PlotManager::~PlotManager() {
    // Don't delete fChain here - it's managed externally
}

std::tuple<Double_t, Double_t, Int_t> PlotManager::getEexcBinningFromIni(
    Double_t defaultMin, Double_t defaultStop, Double_t defaultBin) {
    Double_t eMin = defaultMin;
    Double_t eStop = defaultStop;
    Double_t eBin = defaultBin;

    IniParser inip;
    if (inip.loadFile("reac_info.txt")) {
        if (inip.hasKey("recoil_info", "excEn_start")) eMin = inip.getDouble("recoil_info", "excEn_start");
        if (inip.hasKey("recoil_info", "excEn_stop")) eStop = inip.getDouble("recoil_info", "excEn_stop");
        if (inip.hasKey("recoil_info", "excEn_bin")) eBin = inip.getDouble("recoil_info", "excEn_bin");
    }

    if (eBin <= 0.0) eBin = defaultBin;
    Double_t eMax = eStop + eBin;  // include top generated E* bin edge
    Int_t nBins = static_cast<Int_t>(std::round((eMax - eMin) / eBin));
    if (nBins < 1) nBins = 1;
    return std::make_tuple(eMin, eMax, nBins);
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
    for (Long64_t i = 0; i < nEntries; ++i) {
        fChain->GetEntry(i);
        fillEvent(ejectile, residue, decay_channel);
    }
    
    delete ejectile;
    delete residue;
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
