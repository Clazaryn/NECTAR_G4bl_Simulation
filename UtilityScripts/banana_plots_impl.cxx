// Implementation file for banana plots
#include "make_plots.h"
#include <TMath.h>
#include <cmath>

// ========= BananaPlotManager Implementation =========

BananaPlotManager::BananaPlotManager(TChain* chain, const std::string& setup, 
                                     const std::string& reac, const std::string& rec)
    : PlotManager(chain, setup, reac, rec),
      nBins_dE(100), nBins_E1(250), nBins_Eres(1000),
      dE_max(10), E1_max(25), Eres_max(100) {
    
    // Initialize theta bin edges (0-180 degrees in 10 degree bins)
    for (Int_t i = 0; i <= nThetaBins; ++i) {
        thetaBinEdges[i] = i * 10.0;
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
        
        name = Form("hDEvEres_primary_%s", reaction.c_str());
        title = Form("(dE+E1) vs Eres - Primary Telescope (%s); Eres (MeV); dE+E1 (MeV)", 
                     reaction.c_str());
        hDEvEres_primary.push_back(new TH2D(name, title, nBins_Eres, 0, Eres_max, 
                                             nBins_dE+nBins_E1, 0, dE_max+E1_max));
        
        // Auxillary telescope plots
        name = Form("hDEvE1_auxillary_%s", reaction.c_str());
        title = Form("dE vs E1 - Auxillary Telescope (%s); E1 (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvE1_auxillary.push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                            nBins_dE, 0, dE_max));
        
        name = Form("hDEvEres_auxillary_%s", reaction.c_str());
        title = Form("(dE+E1) vs Eres - Auxillary Telescope (%s); Eres (MeV); dE+E1 (MeV)", 
                     reaction.c_str());
        hDEvEres_auxillary.push_back(new TH2D(name, title, nBins_Eres, 0, Eres_max, 
                                              nBins_dE+nBins_E1, 0, dE_max+E1_max));
        
        // Theta-binned plots for primary telescope
        hDEvE1_primary_theta.resize(nThetaBins);
        hDEvEres_primary_theta.resize(nThetaBins);
        for (Int_t i = 0; i < nThetaBins; ++i) {
            name = Form("hDEvE1_primary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges[i], reaction.c_str());
            title = Form("dE vs E1 - Primary Telescope, #theta = %.0f#circ - %.0f#circ (%s); E1 (MeV); dE (MeV)", 
                        thetaBinEdges[i], thetaBinEdges[i+1], reaction.c_str());
            hDEvE1_primary_theta[i].push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                                        nBins_dE, 0, dE_max));
            
            name = Form("hDEvEres_primary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges[i], reaction.c_str());
            title = Form("(dE+E1) vs Eres - Primary Telescope, #theta = %.0f#circ - %.0f#circ (%s); Eres (MeV); dE+E1 (MeV)", 
                        thetaBinEdges[i], thetaBinEdges[i+1], reaction.c_str());
            hDEvEres_primary_theta[i].push_back(new TH2D(name, title, nBins_Eres, 0, Eres_max, 
                                                         nBins_dE+nBins_E1, 0, dE_max+E1_max));
        }
        
        // Theta-binned plots for auxillary telescope
        hDEvE1_auxillary_theta.resize(nThetaBins);
        hDEvEres_auxillary_theta.resize(nThetaBins);
        for (Int_t i = 0; i < nThetaBins; ++i) {
            name = Form("hDEvE1_auxillary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges[i], reaction.c_str());
            title = Form("dE vs E1 - Auxillary Telescope, #theta = %.0f#circ - %.0f#circ (%s); E1 (MeV); dE (MeV)", 
                        thetaBinEdges[i], thetaBinEdges[i+1], reaction.c_str());
            hDEvE1_auxillary_theta[i].push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                                         nBins_dE, 0, dE_max));
            
            name = Form("hDEvEres_auxillary_theta_%02d_%s", 
                       (Int_t)thetaBinEdges[i], reaction.c_str());
            title = Form("(dE+E1) vs Eres - Auxillary Telescope, #theta = %.0f#circ - %.0f#circ (%s); Eres (MeV); dE+E1 (MeV)", 
                        thetaBinEdges[i], thetaBinEdges[i+1], reaction.c_str());
            hDEvEres_auxillary_theta[i].push_back(new TH2D(name, title, nBins_Eres, 0, Eres_max, 
                                                          nBins_dE+nBins_E1, 0, dE_max+E1_max));
        }
        
    } else if (det_setup == "PoP" || det_setup == "pop") {
        // PoP setup: single telescope
        name = Form("hDEvE1_PoP_%s", reaction.c_str());
        title = Form("dE vs E1 - PoP Telescope (%s); E1 (MeV); dE (MeV)", 
                     reaction.c_str());
        hDEvE1_PoP.push_back(new TH2D(name, title, nBins_E1, 0, E1_max, 
                                      nBins_dE, 0, dE_max));
        
        name = Form("hDEvEres_PoP_%s", reaction.c_str());
        title = Form("(dE+E1) vs Eres - PoP Telescope (%s); Eres (E2+E3+...) (MeV); dE+E1 (MeV)", 
                     reaction.c_str());
        hDEvEres_PoP.push_back(new TH2D(name, title, nBins_Eres, 0, Eres_max, 
                                        nBins_dE+nBins_E1, 0, dE_max+E1_max));
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
    
    // Calculate Eres based on detector setup
    if (det_setup == "new" || det_setup == "New") {
        // New setup: Eres is the Eres detector reading
        Eres = ejectile->meas_Eres;
    } else if (det_setup == "PoP" || det_setup == "pop") {
        // PoP setup: sum E2 + E3 + E4 + E5 + E6
        Eres = ejectile->meas_E2 + ejectile->meas_E3 + ejectile->meas_E4 + ejectile->meas_E5 + ejectile->meas_E6;
    } else {
        Eres = ejectile->meas_Eres;  // Fallback
    }
    
    // Suppress unused parameter warnings
    (void)residue;
    (void)decay_channel;
    
    if (det_setup == "new" || det_setup == "New") {
        if (detector_id == 0) {
            // Primary telescope
            if (hDEvE1_primary.size() > 0) {
                hDEvE1_primary[0]->Fill(E1, dE);
            }
            if (hDEvEres_primary.size() > 0) {
                // For New: (dE+E1) vs Eres (Eres is the Eres detector reading)
                hDEvEres_primary[0]->Fill(Eres, dE + E1);
            }
            
            // Fill theta-binned plots
            Int_t thetaBin = getThetaBin(recon_theta);
            if (thetaBin >= 0 && thetaBin < nThetaBins) {
                size_t thetaBin_size = static_cast<size_t>(thetaBin);
                if (hDEvE1_primary_theta.size() > thetaBin_size && 
                    hDEvE1_primary_theta[thetaBin_size].size() > 0) {
                    hDEvE1_primary_theta[thetaBin_size][0]->Fill(E1, dE);
                }
                if (hDEvEres_primary_theta.size() > thetaBin_size && 
                    hDEvEres_primary_theta[thetaBin_size].size() > 0) {
                    // For New: (dE+E1) vs Eres
                    hDEvEres_primary_theta[thetaBin_size][0]->Fill(Eres, dE + E1);
                }
            }
        } else if (detector_id == 1) {
            // Auxillary telescope
            if (hDEvE1_auxillary.size() > 0) {
                hDEvE1_auxillary[0]->Fill(E1, dE);
            }
            if (hDEvEres_auxillary.size() > 0) {
                // For New: (dE+E1) vs Eres
                hDEvEres_auxillary[0]->Fill(Eres, dE + E1);
            }
            
            // Fill theta-binned plots
            Int_t thetaBin = getThetaBin(recon_theta);
            if (thetaBin >= 0 && thetaBin < nThetaBins) {
                size_t thetaBin_size = static_cast<size_t>(thetaBin);
                if (hDEvE1_auxillary_theta.size() > thetaBin_size && 
                    hDEvE1_auxillary_theta[thetaBin_size].size() > 0) {
                    hDEvE1_auxillary_theta[thetaBin_size][0]->Fill(E1, dE);
                }
                if (hDEvEres_auxillary_theta.size() > thetaBin_size && 
                    hDEvEres_auxillary_theta[thetaBin_size].size() > 0) {
                    // For New: (dE+E1) vs Eres
                    hDEvEres_auxillary_theta[thetaBin_size][0]->Fill(Eres, dE + E1);
                }
            }
        }
    } else if (det_setup == "PoP" || det_setup == "pop") {
        if (hDEvE1_PoP.size() > 0) {
            hDEvE1_PoP[0]->Fill(E1, dE);
        }
        if (hDEvEres_PoP.size() > 0) {
            hDEvEres_PoP[0]->Fill(Eres, dE + E1);
        }
    }
}

Int_t BananaPlotManager::getThetaBin(Double_t theta) const {
    if (theta < thetaBinEdges[0] || theta >= thetaBinEdges[nThetaBins]) {
        return -1;  // Out of range
    }
    
    for (Int_t i = 0; i < nThetaBins; ++i) {
        if (theta >= thetaBinEdges[i] && theta < thetaBinEdges[i+1]) {
            return i;
        }
    }
    
    return -1;  // Should not reach here
}

void BananaPlotManager::writePlots(const char* reaction_name) {
    // Create single output file for all channels combined
    TString output_filename = Form("../%s_sim/plots_%s.root", reaction_name, reaction_name);
    TFile* output_file = new TFile(output_filename.Data(), "RECREATE");
    
    if (!output_file || !output_file->IsOpen()) {
        std::cerr << "banana_plots_impl.cxx:" << __LINE__ << ": Error: Could not create output file: " << output_filename << std::endl;
        return;
    }
    
    output_file->cd();
    
    if (det_setup == "new" || det_setup == "New") {
        // Write primary telescope plots
        for (auto* hist : hDEvE1_primary) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvEres_primary) {
            if (hist) hist->Write();
        }
        
        // Write auxillary telescope plots
        for (auto* hist : hDEvE1_auxillary) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvEres_auxillary) {
            if (hist) hist->Write();
        }
        
        // Create directory for theta-binned plots
        TDirectory* thetaDir = output_file->mkdir("theta_binned");
        thetaDir->cd();
        
        // Write theta-binned plots
        for (Int_t i = 0; i < nThetaBins; ++i) {
            for (auto* hist : hDEvE1_primary_theta[i]) {
                if (hist) hist->Write();
            }
            for (auto* hist : hDEvEres_primary_theta[i]) {
                if (hist) hist->Write();
            }
            for (auto* hist : hDEvE1_auxillary_theta[i]) {
                if (hist) hist->Write();
            }
            for (auto* hist : hDEvEres_auxillary_theta[i]) {
                if (hist) hist->Write();
            }
        }
        
    } else if (det_setup == "PoP" || det_setup == "pop") {
        // Write PoP telescope plots
        for (auto* hist : hDEvE1_PoP) {
            if (hist) hist->Write();
        }
        for (auto* hist : hDEvEres_PoP) {
            if (hist) hist->Write();
        }
    }
    
    output_file->Close();
    delete output_file;
    
    std::cout << "Plots written to: " << output_filename << std::endl;
}
