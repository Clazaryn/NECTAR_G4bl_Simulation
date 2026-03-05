// Implementation file for transmission plots (HR coincidence with virtual detectors).
// Coincidence 2D scatter: MagSept, HRplane, QuadWall. Transmission fraction vs E*.
// New: primary/auxillary. PoP: single. Channels: gamma(red), 1n(blue), 2n(green), 3n(purple), 4n(cyan).
#include "make_plots.h"
#include "ini_parser.h"
#include <TBox.h>
#include <TEllipse.h>
#include <TGraphErrors.h>
#include <set>
#include <sstream>
#include <TLegendEntry.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>

static const Int_t nChannels = 5;
static const Double_t s_coincBinWidth = 2.0;  // 2 MeV bins for coincidence plots
static const Int_t s_maxCoincBins = 6;       // max 2 MeV bins in any channel range
static const char* channelLabels[] = {"#gamma", "1n", "2n", "3n", "4n"};
static const char* channelLabelsROOT[] = {"HRg", "HR1n", "HR2n", "HR3n", "HR4n"};

// Base colors per channel; colour gradient spans channel range (first bin = offset[0])
static const Int_t channelBaseColor[] = {kRed, kBlue, kGreen, kMagenta, kCyan};
static const Int_t channelColorOffsets[s_maxCoincBins] = {3, 2, 1, 0, -4, -7};  // dark to light across 2 MeV bins within each channel's range

// Energy range from reac_info (set in initializePlots)
static Double_t s_Eexc_min = 0.0;
static Double_t s_Eexc_max = 30.0;

// Per-channel energy ranges from calc_hr_ranges logic (set in initializePlots)
static Double_t s_channelStartE[nChannels];
static Double_t s_channelStopE[nChannels];
static Int_t s_channelNBins[nChannels];
static bool s_channelActive[nChannels];

// HR detector box (mm): x in [-137, -15], y in [-20, 20] — matches hrBox in writePlots
static const Double_t s_HRbox_xlo = -132.0, s_HRbox_xhi = -10.0;
static const Double_t s_HRbox_ylo = -20.0, s_HRbox_yhi = 20.0;
static bool isInHRdetectorBox(Double_t x, Double_t y) {
    return (x >= s_HRbox_xlo && x <= s_HRbox_xhi && y >= s_HRbox_ylo && y <= s_HRbox_yhi);
}

// Get coincidence bin index for (channel, Eexc). Returns 0..nBins-1 or -1 if outside channel range.
static Int_t getChannelCoincBin(Int_t channel, Double_t Eexc) {
    if (channel < 0 || channel >= nChannels || !s_channelActive[channel]) return -1;
    if (Eexc < s_channelStartE[channel] || Eexc > s_channelStopE[channel]) return -1;
    Int_t bin = static_cast<Int_t>((Eexc - s_channelStartE[channel]) / s_coincBinWidth);
    if (bin >= s_channelNBins[channel]) bin = s_channelNBins[channel] - 1;
    if (bin < 0) return -1;
    return bin;
}

static Int_t getChannelColor(Int_t channel, Int_t binIndex) {
    if (binIndex < 0) return channelBaseColor[channel];
    Int_t base = channelBaseColor[channel];
    Int_t offIdx = (binIndex >= s_maxCoincBins) ? (s_maxCoincBins - 1) : binIndex;  // clamp to highest colour if > 6
    Int_t off = channelColorOffsets[offIdx];
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

// Compute channel energy ranges (mirrors calc_hr_ranges.py logic)
static void initChannelRanges(const ReactionInfo& rinfo, const std::set<Int_t>& activeChannels) {
    const Double_t overlap_range = 3.0;
    Double_t Sn_CN = rinfo.Sn_CN;
    Double_t Sn_1nDght = rinfo.Sn_1nDght;
    Double_t Sn_2nDght = rinfo.has_Sn_2nDght ? rinfo.Sn_2nDght : 0.0;
    Double_t Sn_3nDght = rinfo.has_Sn_3nDght ? rinfo.Sn_3nDght : 0.0;
    Double_t Sn_4nDght = 0.0;
    IniParser inip;
    if (inip.loadFile("reac_info.txt") && inip.hasKey("separation_energies", "Sn_4nDght"))
        Sn_4nDght = inip.getDouble("separation_energies", "Sn_4nDght");

    for (Int_t ch = 0; ch < nChannels; ++ch) {
        s_channelActive[ch] = (activeChannels.count(ch) > 0);
        s_channelStartE[ch] = s_Eexc_min;
        s_channelStopE[ch] = s_Eexc_min;
    }
    if (activeChannels.count(0)) {
        s_channelStartE[0] = s_Eexc_min;
        s_channelStopE[0] = Sn_CN + overlap_range;
    }
    if (activeChannels.count(1)) {
        s_channelStartE[1] = Sn_CN;
        s_channelStopE[1] = Sn_CN + Sn_1nDght + overlap_range;
    }
    if (activeChannels.count(2)) {
        s_channelStartE[2] = Sn_CN + Sn_1nDght;
        s_channelStopE[2] = Sn_CN + Sn_1nDght + Sn_2nDght + overlap_range;
    }
    if (activeChannels.count(3)) {
        s_channelStartE[3] = Sn_CN + Sn_1nDght + Sn_2nDght;
        s_channelStopE[3] = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + overlap_range;
    }
    if (activeChannels.count(4)) {
        s_channelStartE[4] = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght;
        s_channelStopE[4] = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + Sn_4nDght + overlap_range;
    }
    for (Int_t ch = 0; ch < nChannels; ++ch) {
        s_channelStopE[ch] = std::min(s_channelStopE[ch], s_Eexc_max);
        s_channelNBins[ch] = static_cast<Int_t>(std::ceil((s_channelStopE[ch] - s_channelStartE[ch]) / s_coincBinWidth));
        if (s_channelNBins[ch] < 1) s_channelNBins[ch] = 1;
        if (s_channelNBins[ch] > s_maxCoincBins) s_channelNBins[ch] = s_maxCoincBins;
    }
}

void TransmissionPlotManager::initializePlots() {
    ReactionInfo rinfo;
    rinfo.loadFromIni("reac_info.txt");

    IniParser inip;
    if (inip.loadFile("reac_info.txt")) {
        s_Eexc_min = inip.hasKey("recoil_info", "excEn_start") ? inip.getDouble("recoil_info", "excEn_start") : 0.0;
        s_Eexc_max = inip.hasKey("recoil_info", "excEn_stop") ? inip.getDouble("recoil_info", "excEn_stop") : 30.0;
    }

    // Read excEn_bin for histogram binning (1 MeV default, matching file binning)
    Double_t excEn_bin = 1.0;
    if (inip.hasKey("recoil_info", "excEn_bin"))
        excEn_bin = inip.getDouble("recoil_info", "excEn_bin");
    Int_t nBinsEexc = static_cast<Int_t>(std::round((s_Eexc_max - s_Eexc_min) / excEn_bin));
    if (nBinsEexc < 1) nBinsEexc = 1;

    std::set<Int_t> activeChannels = getActiveChannels();
    if (activeChannels.empty())
        for (Int_t ch = 0; ch < nChannels; ++ch) activeChannels.insert(ch);
    initChannelRanges(rinfo, activeChannels);

    auto createCoincSet = [&](const char* detLabel) {
        std::vector<std::vector<std::vector<TH2D*>>> v(3);
        Double_t xRanges[3][2] = {{-150, 150}, {-150, 0}, {15000, 22000}};  // QuadWall z (mm) absolute G4bl position
        for (Int_t virt = 0; virt < 3; ++virt) {
            v[virt].resize(nChannels);
            const char* vname[] = {"MagSept", "HRplane", "QuadWall"};
            for (Int_t ch = 0; ch < nChannels; ++ch) {
                Int_t nBins = s_channelNBins[ch];
                v[virt][ch].resize(static_cast<size_t>(nBins));
                for (Int_t eb = 0; eb < nBins; ++eb) {
                    Double_t eLo = s_channelStartE[ch] + eb * s_coincBinWidth;
                    Double_t eHi = s_channelStartE[ch] + (eb + 1) * s_coincBinWidth;
                    const char* chLab = (ch < 5) ? channelLabelsROOT[ch] : "ch";
                    TString name = Form("h_%s_%s_%s_eb%d_%s", vname[virt], detLabel, chLab, eb, reaction.c_str());
                    const char* xlabel = (virt == 2) ? "z (mm)" : "x (mm)";  // QuadWall uses (z,y)
                    TString title = Form("%s %s (%.0f-%.0f MeV); %s; y (mm)", vname[virt], channelLabels[ch], eLo, eHi, xlabel);
                    v[virt][ch][static_cast<size_t>(eb)] = new TH2D(name, title, 150, xRanges[virt][0], xRanges[virt][1], 100, -50, 50);
                }
            }
        }
        return v;
    };

    auto createTransmissionSet = [&](const char* detLabel) {
        std::pair<std::vector<TH1D*>, std::vector<TH1D*>> p;
        for (Int_t ch = 0; ch < nChannels; ++ch) {
            const char* chLab = (ch < 5) ? channelLabelsROOT[ch] : "ch";
            TString n1 = Form("h_MagSept_%s_%s_%s", detLabel, chLab, reaction.c_str());
            TString n2 = Form("h_HRplane_%s_%s_%s", detLabel, chLab, reaction.c_str());
            p.first.push_back(new TH1D(n1, Form("MagSept %s; E* (MeV); Counts", channelLabels[ch]), nBinsEexc, s_Eexc_min, s_Eexc_max));
            p.second.push_back(new TH1D(n2, Form("HRplane %s; E* (MeV); Counts", channelLabels[ch]), nBinsEexc, s_Eexc_min, s_Eexc_max));
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
    Int_t coincBin = getChannelCoincBin(decay_channel, Eexc);
    Int_t detId = ejectile->detector_id;

    auto fillCoinc = [&](std::vector<std::vector<std::vector<TH2D*>>>& v) {
        if (coincBin < 0) return;  // Eexc outside channel range
        if (residue->hit_MagSept && v.size() > 0) v[0][decay_channel][static_cast<size_t>(coincBin)]->Fill(residue->MagSept_x, residue->MagSept_y);
        if (residue->hit_HRplane && v.size() > 1) v[1][decay_channel][static_cast<size_t>(coincBin)]->Fill(residue->HRplane_x, residue->HRplane_y);
        if (residue->hit_QuadWall && v.size() > 2) v[2][decay_channel][static_cast<size_t>(coincBin)]->Fill(residue->QuadWall_z, residue->QuadWall_y);
    };

    auto fillTrans = [&](std::vector<TH1D*>& hMag, std::vector<TH1D*>& hHR) {
        if (Eexc < s_Eexc_min || Eexc > s_Eexc_max) return;
        if (residue->hit_MagSept && decay_channel < (Int_t)hMag.size()) hMag[decay_channel]->Fill(Eexc);
        if (residue->hit_HRplane && isInHRdetectorBox(residue->HRplane_x, residue->HRplane_y) && decay_channel < (Int_t)hHR.size())
            hHR[decay_channel]->Fill(Eexc);
    };

    if (det_setup == "new" || det_setup == "New") {
        if (detId == 0) { fillCoinc(h_coinc_primary); fillTrans(h_MagSept_primary, h_HRplane_primary); }
        else if (detId == 1) { fillCoinc(h_coinc_auxillary); fillTrans(h_MagSept_auxillary, h_HRplane_auxillary); }
    } else if (det_setup == "PoP" || det_setup == "pop") {
        fillCoinc(h_coinc_PoP); fillTrans(h_MagSept_PoP, h_HRplane_PoP);
    }
}

// Draw coincidence canvas with boxes (optional beam pipe: single TEllipse arc)
static void drawCoincCanvas(TCanvas*& c, const char* detLabel, std::vector<std::vector<std::vector<TH2D*>>>& h_coinc,
                            Int_t virt, Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi,
                            TBox* box1, TBox* box2, const char* leg1, const char* leg2,
                            const std::set<Int_t>& activeChannels,
                            TEllipse* pipeArc = nullptr, const char* pipeLeg = nullptr) {
    const char* vname[] = {"MagSept", "HRplane", "QuadWall"};
    c = new TCanvas(Form("c_%s_%s", vname[virt], detLabel), Form("%s %s", vname[virt], detLabel), 800, 600);
    c->cd();
    TLegend* leg = new TLegend(0.12, 0.72, 0.22, 0.92);

    TH2D* hFirst = nullptr;
    for (Int_t ch = 0; ch < nChannels; ++ch) {
        if (activeChannels.count(ch) == 0) continue;
        for (Int_t eb = 0; eb < s_channelNBins[ch]; ++eb) {
            TH2D* h = h_coinc[virt][ch][static_cast<size_t>(eb)];
            if (h && h->GetEntries() >= 1) { hFirst = h; break; }
        }
        if (hFirst) break;
    }
    if (!hFirst) { leg->Draw(); c->Update(); return; }

    const char* xlabel = (virt == 2) ? "z (mm)" : "x (mm)";
    hFirst->SetTitle(Form("%s Coincidences (%s); %s; y (mm)", vname[virt], detLabel, xlabel));
    hFirst->GetXaxis()->SetRangeUser(xlo, xhi);
    hFirst->GetYaxis()->SetRangeUser(ylo, yhi);

    // Use a clone only for the axis/frame so the same histogram is never in the pad list twice (avoids TList::Clear "already deleted" on canvas close)
    TH2D* hAxis = static_cast<TH2D*>(hFirst->Clone(Form("hAxis_%s_%s", vname[virt], detLabel)));
    hAxis->SetDirectory(nullptr);
    hAxis->Draw("AXIS");
    if (pipeArc) pipeArc->Draw("SAME");  // after axis so scatter will be on top

    for (Int_t ch = 0; ch < nChannels; ++ch) {
        if (activeChannels.count(ch) == 0) continue;
        bool addedCh = false;
        for (Int_t eb = 0; eb < s_channelNBins[ch]; ++eb) {
            TH2D* h = h_coinc[virt][ch][static_cast<size_t>(eb)];
            if (!h || h->GetEntries() < 1) continue;
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.3);
            h->SetMarkerColor(getChannelColor(ch, eb));
            h->Draw("SCAT SAME");  // each histogram drawn once
            if (!addedCh) {
                TGraph* legDummy = new TGraph(1);
                legDummy->SetPoint(0, 0, 0);
                legDummy->SetMarkerStyle(20);
                legDummy->SetMarkerSize(1.5);
                legDummy->SetMarkerColor(channelBaseColor[ch]);
                leg->AddEntry(legDummy, channelLabels[ch], "p");
                addedCh = true;
            }
        }
    }
    if (box1) { box1->Draw("SAME"); if (leg1) leg->AddEntry(box1, leg1, "l"); }
    if (box2) { box2->Draw("SAME"); if (leg2) leg->AddEntry(box2, leg2, "l"); }
    if (pipeArc && pipeLeg) leg->AddEntry(pipeArc, pipeLeg, "l");
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
    mg->SetTitle(Form("HR transmission fraction (%s); true E* (MeV); HRplane / MagSept (%%)", detLabel));
    mg->GetYaxis()->SetRangeUser(0, 105);
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

        // Beam pipe: full circle r=100 mm, 10 px thick, kGray (one per canvas to avoid TList ownership issues)
        TEllipse* pipeMag = new TEllipse(0, 0, 100, 100);
        pipeMag->SetLineWidth(10); pipeMag->SetLineColor(kGray); pipeMag->SetFillStyle(0);
        TEllipse* pipeHR = new TEllipse(0, 0, 100, 100);
        pipeHR->SetLineWidth(10); pipeHR->SetLineColor(kGray); pipeHR->SetFillStyle(0);

        // Boxes for MagSept (dipole entrance, magnet septum with hashed fill)
        TBox* dipBox = new TBox(-100, -35, 100, 35);
        dipBox->SetLineColor(kRed); dipBox->SetLineWidth(3); dipBox->SetFillStyle(0);
        TBox* sepBox = new TBox(70, -50, 150, 50);
        sepBox->SetLineColor(kBlack); sepBox->SetLineStyle(1); sepBox->SetLineWidth(2);
        sepBox->SetFillStyle(3354); sepBox->SetFillColor(kBlack);  // hashed fill

        drawCoincCanvas(cMag, folder, h_coinc, 0, -150, 150, -50, 50, dipBox, sepBox, "dip entr", "mag sept", activeChannels,
                        pipeMag, "beam pipe");
        if (cMag) cMag->Write();

        // Boxes for HRplane (HR detector position)
        TBox* hrBox = new TBox(s_HRbox_xlo, s_HRbox_ylo, s_HRbox_xhi, s_HRbox_yhi);
        hrBox->SetLineColor(kRed); hrBox->SetLineWidth(3); hrBox->SetFillStyle(0);
        drawCoincCanvas(cHR, folder, h_coinc, 1, -150, 0, -30, 30, hrBox, nullptr, "HR detector", nullptr, activeChannels,
                        pipeHR, "beam pipe");
        if (cHR) cHR->Write();

        drawCoincCanvas(cQW, folder, h_coinc, 2, 15000, 22000, -30, 30, nullptr, nullptr, nullptr, nullptr, activeChannels);
        if (cQW) cQW->Write();

        drawTransmissionCanvas(cTrans, folder, hMag, hHR, activeChannels);
        if (cTrans) cTrans->Write();

        // Write intermediate histograms to detector-labelled histograms subfolder
        TDirectory* histDir = sub->mkdir("histograms");
        histDir->cd();
        for (auto& vv : h_coinc)
            for (auto& vvv : vv)
                for (auto* h : vvv)
                    if (h) h->Write();
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
