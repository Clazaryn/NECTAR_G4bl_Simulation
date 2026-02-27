// Implementation file for excitation energy resolution plots.
// New setup: separate primary and auxillary (routed by detector_id).
// PoP setup: single telescope.
// Output: 2D resolution maps, TGraphErrors sigma vs E*, and TCanvas with comparison plot in exc_resolution/.
#include "make_plots.h"
#include <cmath>
#include <vector>

static const Int_t nBins_Eexc = 26;
static const Double_t Eexc_min = 0.0, Eexc_max = 26.0;
static const Int_t nBins_dEexc = 100;
static const Double_t dEexc_min = -5.0, dEexc_max = 5.0;
static const Int_t minEntriesForFit = 10;

// Helper: fit Gaussian per E* slice and return TGraphErrors
static TGraphErrors* fitResolutionGraph(TH2D* h2d, const char* name, const char* title) {
    std::vector<double> energies, sigmas, energy_errors, sigma_errors;

    for (Int_t xbin = 1; xbin <= h2d->GetNbinsX(); ++xbin) {
        Double_t Eexc_center = h2d->GetXaxis()->GetBinCenter(xbin);
        TH1D* proj = h2d->ProjectionY("_py", xbin, xbin);

        if (!proj || proj->GetEntries() < minEntriesForFit) {
            if (proj) delete proj;
            continue;
        }

        TF1* gausFit = new TF1("gausFit", "gaus", dEexc_min, dEexc_max);
        gausFit->SetParameters(proj->GetMaximum(), 0.0, 0.5);
        Int_t fitStatus = proj->Fit(gausFit, "RQ0");

        if (fitStatus == 0) {
            Double_t sigma = gausFit->GetParameter(2);
            Double_t dsigma = gausFit->GetParError(2);
            energies.push_back(Eexc_center);
            energy_errors.push_back(0.0);
            sigmas.push_back(sigma);
            sigma_errors.push_back(dsigma);
        }

        delete gausFit;
        delete proj;
    }

    if (energies.empty()) return nullptr;

    TGraphErrors* g = new TGraphErrors(static_cast<Int_t>(energies.size()),
                                       energies.data(), sigmas.data(),
                                       energy_errors.data(), sigma_errors.data());
    g->SetName(name);
    g->SetTitle(title);
    return g;
}

// ========= ExcResolutionPlotManager Implementation =========

ExcResolutionPlotManager::ExcResolutionPlotManager(TChain* chain, const std::string& setup,
                                                     const std::string& reac, const std::string& rec)
    : PlotManager(chain, setup, reac, rec),
      h_dEexc_vs_Eexc_primary(nullptr), h_dEexc_vs_Eexc_auxillary(nullptr), h_dEexc_vs_Eexc_PoP(nullptr),
      g_resolution_primary(nullptr), g_resolution_auxillary(nullptr), g_resolution_PoP(nullptr),
      c_resolution(nullptr) {
}

ExcResolutionPlotManager::~ExcResolutionPlotManager() {
    if (h_dEexc_vs_Eexc_primary) { delete h_dEexc_vs_Eexc_primary; h_dEexc_vs_Eexc_primary = nullptr; }
    if (h_dEexc_vs_Eexc_auxillary) { delete h_dEexc_vs_Eexc_auxillary; h_dEexc_vs_Eexc_auxillary = nullptr; }
    if (h_dEexc_vs_Eexc_PoP) { delete h_dEexc_vs_Eexc_PoP; h_dEexc_vs_Eexc_PoP = nullptr; }
    if (g_resolution_primary) { delete g_resolution_primary; g_resolution_primary = nullptr; }
    if (g_resolution_auxillary) { delete g_resolution_auxillary; g_resolution_auxillary = nullptr; }
    if (g_resolution_PoP) { delete g_resolution_PoP; g_resolution_PoP = nullptr; }
    if (c_resolution) { delete c_resolution; c_resolution = nullptr; }
}

void ExcResolutionPlotManager::initializePlots() {
    TString name, title;
    TString reac(reaction.c_str());

    if (det_setup == "new" || det_setup == "New") {
        name = Form("h_dEexc_vs_Eexc_primary_%s", reaction.c_str());
        title = Form("E* accuracy vs E* - Primary (%s); E^{*}_{true} (MeV); E^{*}_{reco} - E^{*}_{true} (MeV)", reac.Data());
        h_dEexc_vs_Eexc_primary = new TH2D(name, title, nBins_Eexc, Eexc_min, Eexc_max, nBins_dEexc, dEexc_min, dEexc_max);

        name = Form("h_dEexc_vs_Eexc_auxillary_%s", reaction.c_str());
        title = Form("E* accuracy vs E* - Auxillary (%s); E^{*}_{true} (MeV); E^{*}_{reco} - E^{*}_{true} (MeV)", reac.Data());
        h_dEexc_vs_Eexc_auxillary = new TH2D(name, title, nBins_Eexc, Eexc_min, Eexc_max, nBins_dEexc, dEexc_min, dEexc_max);
    } else if (det_setup == "PoP" || det_setup == "pop") {
        name = Form("h_dEexc_vs_Eexc_PoP_%s", reaction.c_str());
        title = Form("E* accuracy vs E* (%s); E^{*}_{true} (MeV); E^{*}_{reco} - E^{*}_{true} (MeV)", reac.Data());
        h_dEexc_vs_Eexc_PoP = new TH2D(name, title, nBins_Eexc, Eexc_min, Eexc_max, nBins_dEexc, dEexc_min, dEexc_max);
    }
}

void ExcResolutionPlotManager::fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) {
    (void)residue;
    (void)decay_channel;

    Double_t true_Eexc = ejectile->true_Eexc;
    Double_t recon_Eexc = ejectile->recon_Eexc;

    if (!std::isfinite(recon_Eexc))
        return;

    Double_t dEexc = recon_Eexc - true_Eexc;

    if (det_setup == "new" || det_setup == "New") {
        Int_t detector_id = ejectile->detector_id;
        if (detector_id == 0 && h_dEexc_vs_Eexc_primary)
            h_dEexc_vs_Eexc_primary->Fill(true_Eexc, dEexc);
        else if (detector_id == 1 && h_dEexc_vs_Eexc_auxillary)
            h_dEexc_vs_Eexc_auxillary->Fill(true_Eexc, dEexc);
    } else if (det_setup == "PoP" || det_setup == "pop") {
        if (h_dEexc_vs_Eexc_PoP) h_dEexc_vs_Eexc_PoP->Fill(true_Eexc, dEexc);
    }
}

void ExcResolutionPlotManager::writePlots(const char* reaction_name) {
    TString reac(reaction.c_str());
    TString baseTitle = Form("E* resolution (#sigma) vs E^{*} (%s); E^{*} (MeV); E* res", reac.Data());

    if (det_setup == "new" || det_setup == "New") {
        if (h_dEexc_vs_Eexc_primary)
            g_resolution_primary = fitResolutionGraph(h_dEexc_vs_Eexc_primary, Form("g_Eexc_resolution_primary_%s", reaction_name), baseTitle.Data());
        if (h_dEexc_vs_Eexc_auxillary)
            g_resolution_auxillary = fitResolutionGraph(h_dEexc_vs_Eexc_auxillary, Form("g_Eexc_resolution_auxillary_%s", reaction_name), baseTitle.Data());

        c_resolution = new TCanvas("canv_Eexc_resolution", "E* resolution - Primary vs Auxillary", 800, 600);
        c_resolution->cd();

        bool hasPrimary = g_resolution_primary && g_resolution_primary->GetN() > 0;
        bool hasAuxillary = g_resolution_auxillary && g_resolution_auxillary->GetN() > 0;

        if (hasPrimary) {
            g_resolution_primary->SetMarkerStyle(20);   // filled circle
            g_resolution_primary->SetMarkerColor(kBlack);
            g_resolution_primary->SetMarkerSize(1.2);
            g_resolution_primary->SetMinimum(0);
            g_resolution_primary->Draw("AP");
        }
        if (hasAuxillary) {
            g_resolution_auxillary->SetMarkerStyle(24);  // open circle
            g_resolution_auxillary->SetMarkerColor(kBlack);
            g_resolution_auxillary->SetMarkerSize(1.2);
            g_resolution_auxillary->SetMinimum(0);
            if (hasPrimary)
                g_resolution_auxillary->Draw("P SAME");
            else
                g_resolution_auxillary->Draw("AP");
        }

        if (hasPrimary || hasAuxillary) {
            TLegend* leg = new TLegend(0.65, 0.75, 0.9, 0.9);
            if (hasPrimary) leg->AddEntry(g_resolution_primary, "Primary", "p");
            if (hasAuxillary) leg->AddEntry(g_resolution_auxillary, "Auxillary", "p");
            leg->Draw();
        }

    } else if (det_setup == "PoP" || det_setup == "pop") {
        if (h_dEexc_vs_Eexc_PoP)
            g_resolution_PoP = fitResolutionGraph(h_dEexc_vs_Eexc_PoP, Form("g_Eexc_resolution_PoP_%s", reaction_name), baseTitle.Data());

        c_resolution = new TCanvas("canv_Eexc_resolution", "E* resolution (PoP)", 800, 600);
        c_resolution->cd();
        if (g_resolution_PoP && g_resolution_PoP->GetN() > 0) {
            g_resolution_PoP->SetMarkerStyle(20);
            g_resolution_PoP->SetMinimum(0);
            g_resolution_PoP->Draw("AP");
        }
    }

    // Write to file
    TString output_filename = Form("./%s_results/plots_%s.root", reaction_name, reaction_name);
    TFile* output_file = TFile::Open(output_filename.Data(), "UPDATE");
    if (!output_file || !output_file->IsOpen()) {
        output_file = TFile::Open(output_filename.Data(), "RECREATE");
    }
    if (!output_file || !output_file->IsOpen()) {
        std::cerr << "exc_resolution_plots_impl.cxx:" << __LINE__ << ": Error: Could not open output file: " << output_filename << std::endl;
        return;
    }

    TDirectory* resDir = output_file->GetDirectory("exc_resolution");
    if (!resDir) resDir = output_file->mkdir("exc_resolution");
    if (resDir) {
        resDir->cd();
        if (h_dEexc_vs_Eexc_primary) h_dEexc_vs_Eexc_primary->Write();
        if (h_dEexc_vs_Eexc_auxillary) h_dEexc_vs_Eexc_auxillary->Write();
        if (h_dEexc_vs_Eexc_PoP) h_dEexc_vs_Eexc_PoP->Write();
        if (g_resolution_primary) g_resolution_primary->Write();
        if (g_resolution_auxillary) g_resolution_auxillary->Write();
        if (g_resolution_PoP) g_resolution_PoP->Write();
        if (c_resolution) c_resolution->Write();
    }

    output_file->Close();
    delete output_file;

    std::cout << "Excitation energy resolution plots written to: " << output_filename << "/exc_resolution" << std::endl;
}
