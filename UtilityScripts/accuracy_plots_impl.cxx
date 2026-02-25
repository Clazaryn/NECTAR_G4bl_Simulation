// Implementation file for accuracy plots (reconstructed vs true ejectile properties).
// New setup: primary and auxillary detectors (routing by detector_id from G4bl).
// PoP setup: single detector. Theta axis ranges below are for the vs-#theta 2D histograms only.
#include "make_plots.h"
#include <TMath.h>
#include <cmath>

// Shared bin ranges: E* 0-26 MeV, Eejc 0-150 MeV, accuracy ±5
static const Int_t nBins_Eexc = 26;
static const Double_t Eexc_min = 0.0, Eexc_max = 26.0;
static const Int_t nBins_Eejc = 300;
static const Double_t Eejc_min = 0.0, Eejc_max = 150.0;
static const Int_t nBins_accuracy = 200;
static const Double_t accuracy_min = -5.0, accuracy_max = 5.0;

// Theta axis ranges for the vs-#theta 2D histograms (primary 35-85, auxillary 5-25, PoP 51-69 deg)
static const Double_t Theta_primary_lo = 35.0, Theta_primary_hi = 85.0;
static const Int_t nBins_Theta_primary = 50;
static const Double_t Theta_auxillary_lo = 5.0,  Theta_auxillary_hi = 25.0;
static const Int_t nBins_Theta_auxillary = 20;
static const Double_t Theta_PoP_lo = 51.0, Theta_PoP_hi = 69.0;
static const Int_t nBins_Theta_PoP = 18;

// ========= AccuracyPlotManager Implementation =========

AccuracyPlotManager::AccuracyPlotManager(TChain* chain, const std::string& setup,
                                         const std::string& reac, const std::string& rec)
    : PlotManager(chain, setup, reac, rec) {
}

AccuracyPlotManager::~AccuracyPlotManager() {
    for (auto* h : h_1d_primary) { if (h) delete h; }
    for (auto* h : h_2d_primary) { if (h) delete h; }
    for (auto* h : h_1d_auxillary) { if (h) delete h; }
    for (auto* h : h_2d_auxillary) { if (h) delete h; }
    for (auto* h : h_1d_PoP) { if (h) delete h; }
    for (auto* h : h_2d_PoP) { if (h) delete h; }
}

// Helper: create and push 3x 1D + 9x 2D histograms for one detector
static void createDetectorHists(const std::string& reaction, const std::string& detLabel,
                                Double_t theta_lo, Double_t theta_hi, Int_t nBins_theta,
                                std::vector<TH1D*>& out_1d, std::vector<TH2D*>& out_2d) {
    TString reac(reaction.c_str());
    TString det(detLabel.c_str());
    TString name, title;

    out_1d.clear();
    out_2d.clear();

    name = Form("h_dEexc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E* accuracy %s (%s); E^{*}_{reco} - E^{*}_{true} (MeV); Counts", det.Data(), reac.Data());
    out_1d.push_back(new TH1D(name, title, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dEejc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E_{ejc} accuracy %s (%s); E_{ejc,reco} - E_{ejc,true} (MeV); Counts", det.Data(), reac.Data());
    out_1d.push_back(new TH1D(name, title, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dTheta_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("#theta accuracy %s (%s); #theta_{reco} - #theta_{true} (deg); Counts", det.Data(), reac.Data());
    out_1d.push_back(new TH1D(name, title, nBins_accuracy, accuracy_min, accuracy_max));

    // 2D: vs Eexc (3)
    name = Form("h_dEexc_vs_Eexc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E* accuracy vs E* %s (%s); E^{*}_{true} (MeV); E^{*}_{reco} - E^{*}_{true} (MeV)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_Eexc, Eexc_min, Eexc_max, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dEejc_vs_Eexc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E_{ejc} accuracy vs E* %s (%s); E^{*}_{true} (MeV); E_{ejc,reco} - E_{ejc,true} (MeV)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_Eexc, Eexc_min, Eexc_max, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dTheta_vs_Eexc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("#theta accuracy vs E* %s (%s); E^{*}_{true} (MeV); #theta_{reco} - #theta_{true} (deg)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_Eexc, Eexc_min, Eexc_max, nBins_accuracy, accuracy_min, accuracy_max));
    // 2D: vs Eejc (3)
    name = Form("h_dEexc_vs_Eejc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E* accuracy vs E_{ejc} %s (%s); E_{ejc,true} (MeV); E^{*}_{reco} - E^{*}_{true} (MeV)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_Eejc, Eejc_min, Eejc_max, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dEejc_vs_Eejc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E_{ejc} accuracy vs E_{ejc} %s (%s); E_{ejc,true} (MeV); E_{ejc,reco} - E_{ejc,true} (MeV)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_Eejc, Eejc_min, Eejc_max, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dTheta_vs_Eejc_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("#theta accuracy vs E_{ejc} %s (%s); E_{ejc,true} (MeV); #theta_{reco} - #theta_{true} (deg)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_Eejc, Eejc_min, Eejc_max, nBins_accuracy, accuracy_min, accuracy_max));
    // 2D: vs Theta (3) — detector-specific theta range
    name = Form("h_dEexc_vs_Theta_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E* accuracy vs #theta %s (%s); #theta_{true} (deg); E^{*}_{reco} - E^{*}_{true} (MeV)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_theta, theta_lo, theta_hi, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dEejc_vs_Theta_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("E_{ejc} accuracy vs #theta %s (%s); #theta_{true} (deg); E_{ejc,reco} - E_{ejc,true} (MeV)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_theta, theta_lo, theta_hi, nBins_accuracy, accuracy_min, accuracy_max));
    name = Form("h_dTheta_vs_Theta_%s_%s", reaction.c_str(), detLabel.c_str());
    title = Form("#theta accuracy vs #theta %s (%s); #theta_{true} (deg); #theta_{reco} - #theta_{true} (deg)", det.Data(), reac.Data());
    out_2d.push_back(new TH2D(name, title, nBins_theta, theta_lo, theta_hi, nBins_accuracy, accuracy_min, accuracy_max));
}

void AccuracyPlotManager::initializePlots() {
    if (det_setup == "new" || det_setup == "New") {
        createDetectorHists(reaction, "primary", Theta_primary_lo, Theta_primary_hi, nBins_Theta_primary,
                            h_1d_primary, h_2d_primary);
        createDetectorHists(reaction, "auxillary", Theta_auxillary_lo, Theta_auxillary_hi, nBins_Theta_auxillary,
                            h_1d_auxillary, h_2d_auxillary);
    } else if (det_setup == "PoP" || det_setup == "pop") {
        createDetectorHists(reaction, "PoP", Theta_PoP_lo, Theta_PoP_hi, nBins_Theta_PoP,
                            h_1d_PoP, h_2d_PoP);
    }
}

static void fillOneDetector(Double_t true_Eexc, Double_t true_Eejc, Double_t true_theta,
                           Double_t dEexc, Double_t dEejc, Double_t dTheta,
                           std::vector<TH1D*>& h1, std::vector<TH2D*>& h2) {
    if (h1.size() != 3 || h2.size() != 9) return;
    h1[0]->Fill(dEexc);  h1[1]->Fill(dEejc);  h1[2]->Fill(dTheta);
    h2[0]->Fill(true_Eexc, dEexc);  h2[1]->Fill(true_Eexc, dEejc);  h2[2]->Fill(true_Eexc, dTheta);
    h2[3]->Fill(true_Eejc, dEexc);  h2[4]->Fill(true_Eejc, dEejc);  h2[5]->Fill(true_Eejc, dTheta);
    h2[6]->Fill(true_theta, dEexc); h2[7]->Fill(true_theta, dEejc); h2[8]->Fill(true_theta, dTheta);
}

void AccuracyPlotManager::fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) {
    (void)residue;
    (void)decay_channel;

    Double_t true_Eexc = ejectile->true_Eexc;
    Double_t true_Eejc = ejectile->true_Eejc;
    Double_t true_theta = ejectile->true_theta;
    Double_t recon_Eexc = ejectile->recon_Eexc;
    Double_t recon_Eejc = ejectile->recon_Eejc;
    Double_t recon_theta = ejectile->recon_theta;

    if (!std::isfinite(recon_Eexc) || !std::isfinite(recon_Eejc) || !std::isfinite(recon_theta))
        return;

    Double_t dEexc = recon_Eexc - true_Eexc;
    Double_t dEejc = recon_Eejc - true_Eejc;
    Double_t dTheta = recon_theta - true_theta;

    if (det_setup == "new" || det_setup == "New") {
        Int_t detector_id = ejectile->detector_id;
        if (detector_id == 0)
            fillOneDetector(true_Eexc, true_Eejc, true_theta, dEexc, dEejc, dTheta, h_1d_primary, h_2d_primary);
        else if (detector_id == 1)
            fillOneDetector(true_Eexc, true_Eejc, true_theta, dEexc, dEejc, dTheta, h_1d_auxillary, h_2d_auxillary);
    } else if (det_setup == "PoP" || det_setup == "pop") {
        fillOneDetector(true_Eexc, true_Eejc, true_theta, dEexc, dEejc, dTheta, h_1d_PoP, h_2d_PoP);
    }
}

void AccuracyPlotManager::writePlots(const char* reaction_name) {
    TString output_filename = Form("./%s_results/plots_%s.root", reaction_name, reaction_name);
    TFile* output_file = TFile::Open(output_filename.Data(), "UPDATE");
    if (!output_file || !output_file->IsOpen()) {
        output_file = TFile::Open(output_filename.Data(), "RECREATE");
    }
    if (!output_file || !output_file->IsOpen()) {
        std::cerr << "accuracy_plots_impl.cxx:" << __LINE__ << ": Error: Could not open output file: " << output_filename << std::endl;
        return;
    }

    TDirectory* accDir = output_file->GetDirectory("accuracy_plots");
    if (!accDir) accDir = output_file->mkdir("accuracy_plots");
    if (!accDir) {
        std::cerr << "accuracy_plots_impl.cxx:" << __LINE__ << ": Error: Could not create accuracy_plots directory" << std::endl;
        output_file->Close();
        delete output_file;
        return;
    }

    auto writeDetector = [&](const char* folder, std::vector<TH1D*>& h1, std::vector<TH2D*>& h2) {
        TDirectory* sub = accDir->GetDirectory(folder);
        if (!sub) sub = accDir->mkdir(folder);
        if (!sub) return;
        sub->cd();
        for (auto* h : h1) { if (h) h->Write(); }
        for (auto* h : h2) { if (h) h->Write(); }
    };

    if (det_setup == "new" || det_setup == "New") {
        writeDetector("primary", h_1d_primary, h_2d_primary);
        writeDetector("auxillary", h_1d_auxillary, h_2d_auxillary);
    } else if (det_setup == "PoP" || det_setup == "pop") {
        writeDetector("PoP", h_1d_PoP, h_2d_PoP);
    }

    output_file->Close();
    delete output_file;
    std::cout << "Accuracy plots written to: " << output_filename << "/accuracy_plots" << std::endl;
}
