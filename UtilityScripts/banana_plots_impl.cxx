// Implementation file for banana plots
#include "make_plots.h"
#include <TMath.h>
#include <cmath>

// ========= BananaPlotManager Implementation =========

BananaPlotManager::BananaPlotManager(TChain* chain, const std::string& setup, 
                                     const std::string& reac, const std::string& rec)
    : PlotManager(chain, setup, reac, rec),
      nBins_dE(150), nBins_E1(250), nBins_Eres(250), nBins_Etot(250),
      dE_max(15), E1_max(25), Eres_max(100), Etot_max(140) {
    
    // Primary: 35-45, 45-55, 55-65, 65-75, 75-85 deg
    thetaBinEdges_primary[0] = 35;
    for (Int_t i = 1; i <= nThetaBins_primary; ++i) {
        thetaBinEdges_primary[i] = 35 + i * 10.0;
    }
    // Auxillary: 5-15, 15-25 deg (same structure as primary)
    thetaBinEdges_auxillary[0] = 5;
    for (Int_t i = 1; i <= nThetaBins_auxillary; ++i) {
        thetaBinEdges_auxillary[i] = 5 + i * 10.0;
    }
}

BananaPlotManager::~BananaPlotManager() {
    // ROOT will handle cleanup of histograms when file is closed
}

void BananaPlotManager::initializePlots() {
    TString name, title;
    
    if (det_setup == "new" || det_setup == "New") {
        // New setup: separate histograms for primary and auxillary telescopes
        
        // Primary telescope plots
        name = Form("hDEvE1_primary_%s", reaction.c_str());
        title = Form("dE vs E1 - Primary Telescope (%s); E1 (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvE1_primary.push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                           nBins_dE, 0, dE_max));
        
        name = Form("hDEvEtot_primary_%s", reaction.c_str());
        title = Form("dE vs E_{tot} - Primary Telescope (%s); E_{tot} (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvEtot_primary.push_back(new TH2D(name, title, nBins_Etot, 0, Etot_max, 
                                             nBins_dE, 0, dE_max));
        
        // Auxillary telescope plots
        name = Form("hDEvE1_auxillary_%s", reaction.c_str());
        title = Form("dE vs E1 - Auxillary Telescope (%s); E1 (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvE1_auxillary.push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                            nBins_dE, 0, dE_max));
        
        name = Form("hDEvEtot_auxillary_%s", reaction.c_str());
        title = Form("dE vs E_{tot} - Auxillary Telescope (%s); E_{tot} (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvEtot_auxillary.push_back(new TH2D(name, title, nBins_Etot, 0, Etot_max, 
                                              nBins_dE, 0, dE_max));
        
        // Theta-binned plots for primary telescope (35-45, 45-55, 55-65, 65-75, 75-85 deg)
        hDEvE1_primary_theta.resize(nThetaBins_primary);
        hDEvEtot_primary_theta.resize(nThetaBins_primary);
        for (Int_t i = 0; i < nThetaBins_primary; ++i) {
            name = Form("hDEvE1_primary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges_primary[i], reaction.c_str());
            title = Form("dE vs E1 - Primary Telescope, #theta = %.0f#circ - %.0f#circ (%s); E1 (MeV); dE (MeV)", 
                        thetaBinEdges_primary[i], thetaBinEdges_primary[i+1], reaction.c_str());
            hDEvE1_primary_theta[i].push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                                        nBins_dE, 0, dE_max));
            
            name = Form("hDEvEtot_primary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges_primary[i], reaction.c_str());
            title = Form("dE vs E_{tot} - Primary Telescope, #theta = %.0f#circ - %.0f#circ (%s); E_{tot} (MeV); dE (MeV)", 
                        thetaBinEdges_primary[i], thetaBinEdges_primary[i+1], reaction.c_str());
            hDEvEtot_primary_theta[i].push_back(new TH2D(name, title, nBins_Etot, 0, Etot_max, 
                                                         nBins_dE, 0, dE_max));
        }
        
        // Theta-binned plots for auxillary telescope (5-15, 15-25 deg)
        hDEvE1_auxillary_theta.resize(nThetaBins_auxillary);
        hDEvEtot_auxillary_theta.resize(nThetaBins_auxillary);
        for (Int_t i = 0; i < nThetaBins_auxillary; ++i) {
            name = Form("hDEvE1_auxillary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges_auxillary[i], reaction.c_str());
            title = Form("dE vs E1 - Auxillary Telescope, #theta = %.0f#circ - %.0f#circ (%s); E1 (MeV); dE (MeV)", 
                        thetaBinEdges_auxillary[i], thetaBinEdges_auxillary[i+1], reaction.c_str());
            hDEvE1_auxillary_theta[i].push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                                         nBins_dE, 0, dE_max));
            
            name = Form("hDEvEtot_auxillary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges_auxillary[i], reaction.c_str());
            title = Form("dE vs E_{tot} - Auxillary Telescope, #theta = %.0f#circ - %.0f#circ (%s); E_{tot} (MeV); dE (MeV)", 
                        thetaBinEdges_auxillary[i], thetaBinEdges_auxillary[i+1], reaction.c_str());
            hDEvEtot_auxillary_theta[i].push_back(new TH2D(name, title, nBins_Etot, 0, Etot_max, 
                                                          nBins_dE, 0, dE_max));
        }
        
    } else if (det_setup == "PoP" || det_setup == "pop") {
        // PoP setup: single telescope
        name = Form("hDEvE1_PoP_%s", reaction.c_str());
        title = Form("dE vs E1 - PoP Telescope (%s); E1 (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvE1_PoP.push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                      nBins_dE, 0, dE_max));
        
        name = Form("hDEvEtot_PoP_%s", reaction.c_str());
        title = Form("dE vs E_{tot} - PoP Telescope (%s); E_{tot} (dE+E1+E2+...) (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvEtot_PoP.push_back(new TH2D(name, title, nBins_Etot, 0, Etot_max, 
                                        nBins_dE, 0, dE_max));
    }
}

void BananaPlotManager::fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) {
    // Get measured energies from ejectile object
    // Note: decay_channel is ignored - all channels are combined
    Float_t dE = ejectile->meas_dE;
    Float_t E1 = ejectile->meas_E1;
    Float_t Eres;
    Int_t detector_id = ejectile->detector_id;
    Double_t recon_theta = ejectile->recon_theta;
    
    // Eres: New = Eres detector; PoP = E2+E3+E4+E5+E6
    if (det_setup == "new" || det_setup == "New") {
        Eres = ejectile->meas_Eres;
    } else if (det_setup == "PoP" || det_setup == "pop") {
        Eres = ejectile->meas_E2 + ejectile->meas_E3 + ejectile->meas_E4 + ejectile->meas_E5 + ejectile->meas_E6;
    } else {
        Eres = ejectile->meas_Eres;
    }
    // Punch-through banana: Etot = dE + E1 + Eres (or dE + E1 + E2+... for PoP)
    Float_t Etot = dE + E1 + Eres;
    // Only fill dE vs Etot when ion punches through dE and records a signal in E1
    bool punch_through_dE_E1 = (dE > 0 && E1 > 0);
    
    // Suppress unused parameter warnings
    (void)residue;
    (void)decay_channel;
    
    if (det_setup == "new" || det_setup == "New") {
        if (detector_id == 0) {
            // Primary telescope
            if (hDEvE1_primary.size() > 0) {
                hDEvE1_primary[0]->Fill(E1, dE);
            }
            if (punch_through_dE_E1 && hDEvEtot_primary.size() > 0) {
                hDEvEtot_primary[0]->Fill(Etot, dE);
            }
            
            // Fill theta-binned plots (primary: 35-45,...,75-85 deg)
            Int_t thetaBin = getThetaBinPrimary(recon_theta);
            if (thetaBin >= 0 && thetaBin < nThetaBins_primary) {
                size_t thetaBin_size = static_cast<size_t>(thetaBin);
                if (hDEvE1_primary_theta.size() > thetaBin_size && 
                    hDEvE1_primary_theta[thetaBin_size].size() > 0) {
                    hDEvE1_primary_theta[thetaBin_size][0]->Fill(E1, dE);
                }
                if (punch_through_dE_E1 && hDEvEtot_primary_theta.size() > thetaBin_size && 
                    hDEvEtot_primary_theta[thetaBin_size].size() > 0) {
                    hDEvEtot_primary_theta[thetaBin_size][0]->Fill(Etot, dE);
                }
            }
        } else if (detector_id == 1) {
            // Auxillary telescope
            if (hDEvE1_auxillary.size() > 0) {
                hDEvE1_auxillary[0]->Fill(E1, dE);
            }
            if (punch_through_dE_E1 && hDEvEtot_auxillary.size() > 0) {
                hDEvEtot_auxillary[0]->Fill(Etot, dE);
            }
            
            // Fill theta-binned plots (auxillary: 5-15, 15-25 deg)
            Int_t thetaBin = getThetaBinAuxillary(recon_theta);
            if (thetaBin >= 0 && thetaBin < nThetaBins_auxillary) {
                size_t thetaBin_size = static_cast<size_t>(thetaBin);
                if (hDEvE1_auxillary_theta.size() > thetaBin_size && 
                    hDEvE1_auxillary_theta[thetaBin_size].size() > 0) {
                    hDEvE1_auxillary_theta[thetaBin_size][0]->Fill(E1, dE);
                }
                if (punch_through_dE_E1 && hDEvEtot_auxillary_theta.size() > thetaBin_size && 
                    hDEvEtot_auxillary_theta[thetaBin_size].size() > 0) {
                    hDEvEtot_auxillary_theta[thetaBin_size][0]->Fill(Etot, dE);
                }
            }
        }
    } else if (det_setup == "PoP" || det_setup == "pop") {
        if (hDEvE1_PoP.size() > 0) {
            hDEvE1_PoP[0]->Fill(E1, dE);
        }
        if (punch_through_dE_E1 && hDEvEtot_PoP.size() > 0) {
            hDEvEtot_PoP[0]->Fill(Etot, dE);
        }
    }
}

Int_t BananaPlotManager::getThetaBinPrimary(Double_t theta) const {
    if (theta < thetaBinEdges_primary[0] || theta >= thetaBinEdges_primary[nThetaBins_primary]) {
        return -1;
    }
    for (Int_t i = 0; i < nThetaBins_primary; ++i) {
        if (theta >= thetaBinEdges_primary[i] && theta < thetaBinEdges_primary[i+1]) {
            return i;
        }
    }
    return -1;
}

Int_t BananaPlotManager::getThetaBinAuxillary(Double_t theta) const {
    if (theta < thetaBinEdges_auxillary[0] || theta >= thetaBinEdges_auxillary[nThetaBins_auxillary]) {
        return -1;
    }
    for (Int_t i = 0; i < nThetaBins_auxillary; ++i) {
        if (theta >= thetaBinEdges_auxillary[i] && theta < thetaBinEdges_auxillary[i+1]) {
            return i;
        }
    }
    return -1;
}

void BananaPlotManager::writePlots(const char* reaction_name) {
    TString output_filename = Form("./%s_results/plots_%s.root", reaction_name, reaction_name);
    TFile* output_file = new TFile(output_filename.Data(), "RECREATE");
    
    if (!output_file || !output_file->IsOpen()) {
        std::cerr << "banana_plots_impl.cxx:" << __LINE__ << ": Error: Could not create output file: " << output_filename << std::endl;
        return;
    }
    
    TDirectory* bananaDir = output_file->mkdir("banana_plots");
    if (!bananaDir) {
        std::cerr << "banana_plots_impl.cxx:" << __LINE__ << ": Error: Could not create banana_plots directory" << std::endl;
        output_file->Close();
        delete output_file;
        return;
    }
    bananaDir->cd();
    
    if (det_setup == "new" || det_setup == "New") {
        for (auto* hist : hDEvE1_primary) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvEtot_primary) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvE1_auxillary) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvEtot_auxillary) {
            if (hist) hist->Write();
        }
        
        TDirectory* thetaDir = bananaDir->mkdir("theta_binned");
        if (thetaDir) {
            thetaDir->cd();
            for (Int_t i = 0; i < nThetaBins_primary; ++i) {
                for (auto* hist : hDEvE1_primary_theta[i]) {
                    if (hist) hist->Write();
                }
                for (auto* hist : hDEvEtot_primary_theta[i]) {
                    if (hist) hist->Write();
                }
            }
            for (Int_t i = 0; i < nThetaBins_auxillary; ++i) {
                for (auto* hist : hDEvE1_auxillary_theta[i]) {
                    if (hist) hist->Write();
                }
                for (auto* hist : hDEvEtot_auxillary_theta[i]) {
                    if (hist) hist->Write();
                }
            }
        }
        
    } else if (det_setup == "PoP" || det_setup == "pop") {
        for (auto* hist : hDEvE1_PoP) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvEtot_PoP) {
            if (hist) hist->Write();
        }
    }
    
    output_file->Close();
    delete output_file;
    
    std::cout << "Plots written to: " << output_filename << "/banana_plots" << std::endl;
}
