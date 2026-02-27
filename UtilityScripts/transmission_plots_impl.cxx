// Implementation file for transmission plots (HR coincidence with virtual detectors).
// Coincidence 2D scatter: MagSept, HRplane, QuadWall. Transmission fraction vs E*.
// New: primary/auxillary. PoP: single. Channels: gamma(red), 1n(blue), 2n(green), 3n(purple), 4n(cyan).
#include "make_plots.h"
#include "ini_parser.h"
#include <TBox.h>
#include <TGraphErrors.h>
#include <set>
#include <sstream>
#include <TLegendEntry.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <cmath>
#include <vector>
#include <map>

static const Int_t nChannels = 5;
static const Int_t nEexcBins = 5;
static const Double_t Eexc_edges[nEexcBins + 1] = {0.0, 6.0, 12.0, 18.0, 24.0, 30.0};
static const char* channelLabels[] = {"#gamma", "1n", "2n", "3n", "4n"};

// Base colors per channel (red, blue, green, purple, cyan); hue gradient uses offsets
static const Int_t channelBaseColor[] = {kRed, kBlue, kGreen, kMagenta, kCyan};
static const Int_t channelColorOffsets[] = {-4, -2, 0, 2, 4};  // for Eexc bins: dark to light

static Int_t getEexcBin(Double_t Eexc) {
    for (Int_t i = 0; i < nEexcBins; ++i)
        if (Eexc >= Eexc_edges[i] && Eexc < Eexc_edges[i + 1]) return i;
    if (Eexc >= Eexc_edges[nEexcBins]) return nEexcBins - 1;
    return 0;
}

static Int_t getChannelColor(Int_t channel, Int_t eexcBin) {
    Int_t base = channelBaseColor[channel];
    Int_t off = channelColorOffsets[eexcBin];
    return base + off;
}

static std::set<Int_t> getActiveChannels() {
    std::set<Int_t> active;
    IniParser parser;
    if (!parser.loadFile("reac_info.txt") || !parser.hasKey("recoil_info", "run_HR_modes"))
        return active;
    std::string modes = parser.getString("recoil_info", "run_HR_modes");
    std::istringstream ss(modes);
    std::string mode;
    std::map<std::string, Int_t> modeToCh = {{"HRg", 0}, {"HR1n", 1}, {"HR2n", 2}, {"HR3n", 3}, {"HR4n", 4}};
    while (std::getline(ss, mode, ',')) {
        size_t start = mode.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        size_t end = mode.find_last_not_of(" \t");
        mode = (end == std::string::npos) ? mode.substr(start) : mode.substr(start, end - start + 1);
        auto it = modeToCh.find(mode);
        if (it != modeToCh.end()) active.insert(it->second);
    }
    return active;
}

// ========= TransmissionPlotManager Implementation =========

TransmissionPlotManager::TransmissionPlotManager(TChain* chain, const std::string& setup,
                                                 const std::string& reac, const std::string& rec)
    : PlotManager(chain, setup, reac, rec),
      c_MagSept_primary(nullptr), c_HRplane_primary(nullptr), c_QuadWall_primary(nullptr), c_transmission_primary(nullptr),
      c_MagSept_auxillary(nullptr), c_HRplane_auxillary(nullptr), c_QuadWall_auxillary(nullptr), c_transmission_auxillary(nullptr),
      c_MagSept_PoP(nullptr), c_HRplane_PoP(nullptr), c_QuadWall_PoP(nullptr), c_transmission_PoP(nullptr) {
}

TransmissionPlotManager::~TransmissionPlotManager() {
    auto deleteCoinc = [](std::vector<std::vector<std::vector<TH2D*>>>& v) {
        for (auto& vv : v) for (auto& vvv : vv) for (auto* h : vvv) delete h;
    };
    deleteCoinc(h_coinc_primary); deleteCoinc(h_coinc_auxillary); deleteCoinc(h_coinc_PoP);
    for (auto* h : h_MagSept_primary) delete h; for (auto* h : h_HRplane_primary) delete h;
    for (auto* h : h_MagSept_auxillary) delete h; for (auto* h : h_HRplane_auxillary) delete h;
    for (auto* h : h_MagSept_PoP) delete h; for (auto* h : h_HRplane_PoP) delete h;
    if (c_MagSept_primary) delete c_MagSept_primary; if (c_HRplane_primary) delete c_HRplane_primary; if (c_QuadWall_primary) delete c_QuadWall_primary;
    if (c_MagSept_auxillary) delete c_MagSept_auxillary; if (c_HRplane_auxillary) delete c_HRplane_auxillary; if (c_QuadWall_auxillary) delete c_QuadWall_auxillary;
    if (c_MagSept_PoP) delete c_MagSept_PoP; if (c_HRplane_PoP) delete c_HRplane_PoP; if (c_QuadWall_PoP) delete c_QuadWall_PoP;
    if (c_transmission_primary) delete c_transmission_primary; if (c_transmission_auxillary) delete c_transmission_auxillary; if (c_transmission_PoP) delete c_transmission_PoP;
}

void TransmissionPlotManager::initializePlots() {
    ReactionInfo rinfo;
    rinfo.loadFromIni("reac_info.txt");
    Int_t nEexc = static_cast<Int_t>(rinfo.recoil_excEns.size());
    if (nEexc < 1) nEexc = 31;
    Double_t Eexc_min = rinfo.recoil_excEns.empty() ? 0.0 : rinfo.recoil_excEns.front();
    Double_t Eexc_max = rinfo.recoil_excEns.empty() ? 30.0 : rinfo.recoil_excEns.back();

    auto createCoincSet = [&](const char* detLabel) {
        std::vector<std::vector<std::vector<TH2D*>>> v(3);
        Double_t xRanges[3][2] = {{-150, 150}, {-150, 0}, {2361, 3361}};
        for (Int_t virt = 0; virt < 3; ++virt) {
            v[virt].resize(nChannels);
            const char* vname[] = {"MagSept", "HRplane", "QuadWall"};
            for (Int_t ch = 0; ch < nChannels; ++ch) {
                v[virt][ch].resize(nEexcBins);
                for (Int_t eb = 0; eb < nEexcBins; ++eb) {
                    TString name = Form("h_%s_%s_ch%d_eb%d_%s", vname[virt], detLabel, ch, eb, reaction.c_str());
                    TString title = Form("%s %s (%.0f-%.0f MeV); x (mm); y (mm)", vname[virt], channelLabels[ch], Eexc_edges[eb], Eexc_edges[eb+1]);
                    v[virt][ch][eb] = new TH2D(name, title, 150, xRanges[virt][0], xRanges[virt][1], 100, -50, 50);
                }
            }
        }
        return v;
    };

    auto createTransmissionSet = [&](const char* detLabel) {
        std::pair<std::vector<TH1D*>, std::vector<TH1D*>> p;
        for (Int_t ch = 0; ch < nChannels; ++ch) {
            TString n1 = Form("h_MagSept_%s_ch%d_%s", detLabel, ch, reaction.c_str());
            TString n2 = Form("h_HRplane_%s_ch%d_%s", detLabel, ch, reaction.c_str());
            p.first.push_back(new TH1D(n1, Form("MagSept %s; E* (MeV); Counts", channelLabels[ch]), nEexc, Eexc_min, Eexc_max));
            p.second.push_back(new TH1D(n2, Form("HRplane %s; E* (MeV); Counts", channelLabels[ch]), nEexc, Eexc_min, Eexc_max));
        }
        return p;
    };

    if (det_setup == "new" || det_setup == "New") {
        h_coinc_primary = createCoincSet("primary");
        h_coinc_auxillary = createCoincSet("auxillary");
        auto tp = createTransmissionSet("primary");
        h_MagSept_primary = tp.first; h_HRplane_primary = tp.second;
        auto ta = createTransmissionSet("auxillary");
        h_MagSept_auxillary = ta.first; h_HRplane_auxillary = ta.second;
    } else if (det_setup == "PoP" || det_setup == "pop") {
        h_coinc_PoP = createCoincSet("PoP");
        auto t = createTransmissionSet("PoP");
        h_MagSept_PoP = t.first; h_HRplane_PoP = t.second;
    }
}

void TransmissionPlotManager::fillEvent(LightEjectile* ejectile, HeavyResidue* residue, Short_t decay_channel) {
    if (decay_channel < 0 || decay_channel >= nChannels) return;

    Double_t Eexc = ejectile->true_Eexc;
    Int_t eexcBin = getEexcBin(Eexc);
    Int_t detId = ejectile->detector_id;

    auto fillCoinc = [&](std::vector<std::vector<std::vector<TH2D*>>>& v) {
        if (residue->hit_MagSept && v.size() > 0) v[0][decay_channel][eexcBin]->Fill(residue->MagSept_x, residue->MagSept_y);
        if (residue->hit_HRplane && v.size() > 1) v[1][decay_channel][eexcBin]->Fill(residue->HRplane_x, residue->HRplane_y);
        if (residue->hit_QuadWall && v.size() > 2) v[2][decay_channel][eexcBin]->Fill(residue->QuadWall_x, residue->QuadWall_y);
    };

    auto fillTrans = [&](std::vector<TH1D*>& hMag, std::vector<TH1D*>& hHR) {
        if (residue->hit_MagSept && decay_channel < (Int_t)hMag.size()) hMag[decay_channel]->Fill(Eexc);
        if (residue->hit_HRplane && decay_channel < (Int_t)hHR.size()) hHR[decay_channel]->Fill(Eexc);
    };

    if (det_setup == "new" || det_setup == "New") {
        if (detId == 0) { fillCoinc(h_coinc_primary); fillTrans(h_MagSept_primary, h_HRplane_primary); }
        else if (detId == 1) { fillCoinc(h_coinc_auxillary); fillTrans(h_MagSept_auxillary, h_HRplane_auxillary); }
    } else if (det_setup == "PoP" || det_setup == "pop") {
        fillCoinc(h_coinc_PoP); fillTrans(h_MagSept_PoP, h_HRplane_PoP);
    }
}

// Draw coincidence canvas with boxes
static void drawCoincCanvas(TCanvas*& c, const char* detLabel, std::vector<std::vector<std::vector<TH2D*>>>& h_coinc,
                            Int_t virt, Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi,
                            TBox* box1, TBox* box2, const char* leg1, const char* leg2,
                            const std::set<Int_t>& activeChannels) {
    const char* vname[] = {"MagSept", "HRplane", "QuadWall"};
    c = new TCanvas(Form("c_%s_%s", vname[virt], detLabel), Form("%s %s", vname[virt], detLabel), 800, 600);
    c->cd();
    TLegend* leg = new TLegend(0.78, 0.35, 0.98, 0.65);
    bool first = true;
    for (Int_t ch = 0; ch < nChannels; ++ch) {
        if (activeChannels.count(ch) == 0) continue;
        bool addedCh = false;
        for (Int_t eb = 0; eb < nEexcBins; ++eb) {
            TH2D* h = h_coinc[virt][ch][eb];
            if (!h || h->GetEntries() < 1) continue;
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.3);
            h->SetMarkerColor(getChannelColor(ch, eb));
            if (first) {
                h->Draw("SCAT");
                h->GetXaxis()->SetRangeUser(xlo, xhi);
                h->GetYaxis()->SetRangeUser(ylo, yhi);
                first = false;
            } else {
                h->Draw("SCAT SAME");
            }
            if (!addedCh) {
                TGraph* legDummy = new TGraph(1);
                legDummy->SetPoint(0, 0, 0);
                legDummy->SetMarkerStyle(20);
                legDummy->SetMarkerSize(1.5);
                legDummy->SetMarkerColor(getChannelColor(ch, eb));
                leg->AddEntry(legDummy, channelLabels[ch], "p");
                addedCh = true;
            }
        }
    }
    if (box1) { box1->Draw("SAME"); if (leg1) leg->AddEntry(box1, leg1, "l"); }
    if (box2) { box2->Draw("SAME"); if (leg2) leg->AddEntry(box2, leg2, "l"); }
    leg->Draw();
    c->Update();
}

// Transmission fraction canvas
static void drawTransmissionCanvas(TCanvas*& c, const char* detLabel,
                                   std::vector<TH1D*>& hMag, std::vector<TH1D*>& hHR,
                                   const std::set<Int_t>& activeChannels) {
    c = new TCanvas(Form("c_trans_%s", detLabel), Form("Transmission %s", detLabel), 800, 600);
    c->cd();
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.88, 0.7, 0.98, 0.95);
    for (Int_t ch = 0; ch < nChannels && ch < (Int_t)hMag.size() && ch < (Int_t)hHR.size(); ++ch) {
        if (activeChannels.count(ch) == 0) continue;
        TH1D* hM = hMag[ch];
        TH1D* hH = hHR[ch];
        if (!hM || !hH) continue;
        Int_t n = hM->GetNbinsX();
        std::vector<double> xv, yv, exv, eyv;
        for (Int_t i = 1; i <= n; ++i) {
            Double_t mag = hM->GetBinContent(i);
            Double_t hr = hH->GetBinContent(i);
            if (mag <= 0 || hr <= 0) continue;  // suppress 0% transmission points
            Double_t frac = hr / mag;
            Double_t yval = 100.0 * frac;
            xv.push_back(hM->GetBinCenter(i));
            exv.push_back(0.0);
            yv.push_back(yval);
            eyv.push_back(100.0 * frac * std::sqrt(1.0 / hr + 1.0 / mag));
        }
        if (xv.empty()) continue;
        TGraphErrors* g = new TGraphErrors(static_cast<Int_t>(xv.size()), xv.data(), yv.data(), exv.data(), eyv.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->SetMarkerColor(channelBaseColor[ch]);
        mg->Add(g);
        leg->AddEntry(g, channelLabels[ch], "p");
    }
    mg->Draw("AP");
    mg->SetTitle("HR transmission fraction; E* (MeV); HRplane / MagSept (%)");
    leg->Draw();
    c->Update();
}

void TransmissionPlotManager::writePlots(const char* reaction_name) {
    TString outPath = Form("./%s_results/plots_%s.root", reaction_name, reaction_name);
    TFile* f = TFile::Open(outPath.Data(), "UPDATE");
    if (!f || !f->IsOpen()) f = TFile::Open(outPath.Data(), "RECREATE");
    if (!f || !f->IsOpen()) {
        std::cerr << "transmission_plots_impl.cxx: Could not open " << outPath << std::endl;
        return;
    }

    TDirectory* transDir = f->GetDirectory("transmission_plots");
    if (!transDir) transDir = f->mkdir("transmission_plots");

    std::set<Int_t> activeChannels = getActiveChannels();
    if (activeChannels.empty()) {
        for (Int_t ch = 0; ch < nChannels; ++ch) activeChannels.insert(ch);
    }

    auto writeTelescope = [&](const char* folder, std::vector<std::vector<std::vector<TH2D*>>>& h_coinc,
                             std::vector<TH1D*>& hMag, std::vector<TH1D*>& hHR,
                             TCanvas*& cMag, TCanvas*& cHR, TCanvas*& cQW, TCanvas*& cTrans) {
        TDirectory* sub = transDir->mkdir(folder);
        sub->cd();

        // Boxes for MagSept
        TBox* dipBox = new TBox(-100, -35, 100, 35);
        dipBox->SetLineColor(kRed); dipBox->SetLineWidth(3); dipBox->SetFillStyle(0);
        TBox* sepBox = new TBox(70, -50, 150, 50);
        sepBox->SetLineColor(kBlack); sepBox->SetLineStyle(2); sepBox->SetLineWidth(2); sepBox->SetFillStyle(0);

        drawCoincCanvas(cMag, folder, h_coinc, 0, -150, 150, -50, 50, dipBox, sepBox, "dipole entrance", "magnet septum", activeChannels);
        if (cMag) cMag->Write();

        TBox* hrBox = new TBox(-137, -20, -15, 20);
        hrBox->SetLineColor(kRed); hrBox->SetLineWidth(3); hrBox->SetFillStyle(0);
        drawCoincCanvas(cHR, folder, h_coinc, 1, -150, 0, -30, 30, hrBox, nullptr, "HR detector position", nullptr, activeChannels);
        if (cHR) cHR->Write();

        drawCoincCanvas(cQW, folder, h_coinc, 2, 2361, 3361, -30, 30, nullptr, nullptr, nullptr, nullptr, activeChannels);
        if (cQW) cQW->Write();

        drawTransmissionCanvas(cTrans, folder, hMag, hHR, activeChannels);
        if (cTrans) cTrans->Write();

        for (auto& vv : h_coinc) for (auto& vvv : vv) for (auto* h : vvv) if (h) h->Write();
        for (auto* h : hMag) if (h) h->Write();
        for (auto* h : hHR) if (h) h->Write();
    };

    if (det_setup == "new" || det_setup == "New") {
        writeTelescope("primary", h_coinc_primary, h_MagSept_primary, h_HRplane_primary,
                       c_MagSept_primary, c_HRplane_primary, c_QuadWall_primary, c_transmission_primary);
        writeTelescope("auxillary", h_coinc_auxillary, h_MagSept_auxillary, h_HRplane_auxillary,
                       c_MagSept_auxillary, c_HRplane_auxillary, c_QuadWall_auxillary, c_transmission_auxillary);
        h_coinc_primary.clear(); h_coinc_auxillary.clear();
        h_MagSept_primary.clear(); h_HRplane_primary.clear();
        h_MagSept_auxillary.clear(); h_HRplane_auxillary.clear();
        c_MagSept_primary = c_HRplane_primary = c_QuadWall_primary = c_transmission_primary = nullptr;
        c_MagSept_auxillary = c_HRplane_auxillary = c_QuadWall_auxillary = c_transmission_auxillary = nullptr;
    } else if (det_setup == "PoP" || det_setup == "pop") {
        writeTelescope("PoP", h_coinc_PoP, h_MagSept_PoP, h_HRplane_PoP,
                       c_MagSept_PoP, c_HRplane_PoP, c_QuadWall_PoP, c_transmission_PoP);
        h_coinc_PoP.clear();
        h_MagSept_PoP.clear(); h_HRplane_PoP.clear();
        c_MagSept_PoP = c_HRplane_PoP = c_QuadWall_PoP = c_transmission_PoP = nullptr;
    }

    f->Close();
    delete f;
    std::cout << "Transmission plots written to: " << outPath << "/transmission_plots" << std::endl;
}
