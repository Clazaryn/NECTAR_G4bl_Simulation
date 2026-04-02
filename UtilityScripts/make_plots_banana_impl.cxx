// Implementation file for make_plots_banana.h (base banana class)
#include "make_plots_banana.h"
#include <iostream>

// ========= BananaPlotManagerBase Implementation =========

BananaPlotManagerBase::BananaPlotManagerBase(TChain* chain,
                                             const std::string& setup,
                                             const std::string& reac,
                                             const std::string& rec)
    : fChain(chain), det_setup(setup), reaction(reac), recType(rec) {
}

BananaPlotManagerBase::~BananaPlotManagerBase() {
    // Don't delete fChain here - it's managed externally
}

void BananaPlotManagerBase::fillPlots() {
    if (!fChain) {
        std::cerr << "make_plots_banana_impl.cxx:" << __LINE__
                  << ": Error: fChain is null" << std::endl;
        return;
    }

    // Set up branches
    LightEjectile* ejectile = new LightEjectile();

    fChain->SetBranchAddress("ejectile", &ejectile);

    // Loop over all events
    Long64_t nEntries = fChain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        fChain->GetEntry(i);
        fillEvent(ejectile);
    }

    delete ejectile;
}
