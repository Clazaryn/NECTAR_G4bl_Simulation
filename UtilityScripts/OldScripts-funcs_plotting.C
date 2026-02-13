#include "funcs_plotting.h"

// Utility: Set X and Y axis ranges for first TH2 on a canvas
void SetTH2RangeOnCanvas(TCanvas* canvas, double xmin, double xmax, double ymin, double ymax) {
    if (!canvas) return;
    canvas->cd();
    TList* prims = canvas->GetListOfPrimitives();
    if (!prims) return;
    for (int i = 0; i < prims->GetSize(); ++i) {
        TObject* obj = prims->At(i);
        if (obj && obj->InheritsFrom("TH2")) {
            TH2* h = (TH2*)obj;
            h->GetXaxis()->SetRangeUser(xmin, xmax);
            h->GetYaxis()->SetRangeUser(ymin, ymax);
            canvas->Modified();
            canvas->Update();
            break;
        }
    }
}

// Create a custom color palette from two RGB colors
std::vector<int> CreateCustomPalette(const std::vector<double>& energies,
    double r1, double g1, double b1,  // Start color
    double r2, double g2, double b2   // End color
) {
    int nColors = energies.size();
    std::vector<int> palette;
    std::vector<double> stops(nColors), red(nColors), green(nColors), blue(nColors);

    // Fill gradient arrays
    for (int i = 0; i < nColors; ++i) {
        double t = static_cast<double>(i) / (nColors - 1);
        stops[i] = t;
        red[i]   = (1 - t) * r1 + t * r2;
        green[i] = (1 - t) * g1 + t * g2;
        blue[i]  = (1 - t) * b1 + t * b2;
    }

    int baseIndex = TColor::CreateGradientColorTable(nColors, stops.data(), red.data(), green.data(), blue.data(), nColors);
    for (int i = 0; i < nColors; ++i) { palette.push_back(baseIndex + i); }

    gStyle->SetPalette(nColors, palette.data());
    return palette;
}

// Map excitation energies to color indices
std::map<double, int> MapEnergiesToColors(const std::vector<double>& energies,
    double r1, double g1, double b1,
    double r2, double g2, double b2
) {
    std::vector<double> sortedE = energies;
    std::sort(sortedE.begin(), sortedE.end());

    std::vector<int> colors = CreateCustomPalette(sortedE, r1, g1, b1, r2, g2, b2);
    std::map<double, int> energyColorMap;

    for (size_t i = 0; i < sortedE.size(); ++i) { energyColorMap[sortedE[i]] = colors[i]; }

    return energyColorMap;
}

// Fit Gaussians and generate resolution graphs
TGraphErrors* FitResolutionGraph(TFile* file, const TString& prefix, TH2D* resolutionHist) {
    std::vector<double> energies, fwhms, energy_errors, fwhm_errors;
    std::map<std::string, TKey*> latestKeys;

    TIter next(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) 
    {
        std::string baseName = key->GetName();
        if (latestKeys.find(baseName) == latestKeys.end() || key->GetCycle() > latestKeys[baseName]->GetCycle())
            latestKeys[baseName] = key;
    }

    for (const auto& pair : latestKeys) 
    {
        TKey* key = pair.second;        TObject* obj = key->ReadObj();
        if (!obj->InheritsFrom("TH1D")) continue;

        TH1D* hist = (TH1D*)obj;        TString name = hist->GetName();
        if (!name.BeginsWith(prefix)) continue;     
        if (!name.EndsWith("_Eexc_spect")) continue;
        if (hist->GetEntries() < 10) continue;

        TString nameCopy = name;
        nameCopy.ReplaceAll(prefix, "");    nameCopy.ReplaceAll("_Eexc_spect", "");     nameCopy.ReplaceAll("p", ".");
        double Eexc = nameCopy.Atof();

        TF1* gausFit = new TF1("gausFit", "gaus", Eexc - 2.5, Eexc + 2.5);
        gausFit->SetParameters(50.0, Eexc, 0.5);
        hist->Fit(gausFit, "RQ");

        if (!gausFit) continue;
        double sigma = gausFit->GetParameter(2);    double dsigma = gausFit->GetParError(2);
        energies.push_back(Eexc);       energy_errors.push_back(0);
        fwhms.push_back(sigma);         fwhm_errors.push_back(dsigma);
        
        // Fill TH2D if provided
        if (resolutionHist) {
            // Create a copy of the histogram and subtract the Eexc value
            TH1D* histCopy = (TH1D*)hist->Clone("temp_hist");
            
            // Shift the histogram by subtracting Eexc from all bin centers
            for (int bin = 1; bin <= histCopy->GetNbinsX(); ++bin) {
                double oldCenter = histCopy->GetBinCenter(bin);     double newCenter = oldCenter - Eexc;
                double content = histCopy->GetBinContent(bin);      double error = histCopy->GetBinError(bin);
                
                int newBin = resolutionHist->GetYaxis()->FindBin(newCenter);
                if (newBin > 0 && newBin <= resolutionHist->GetNbinsY()) {
                    int xBin = resolutionHist->GetXaxis()->FindBin(Eexc);
                    if (xBin > 0 && xBin <= resolutionHist->GetNbinsX()) {
                        resolutionHist->SetBinContent(xBin, newBin, content);
                    }
                }
            }
            delete histCopy;
        }
    }
    
    if (energies.empty()) { return nullptr;}
    else { return new TGraphErrors(energies.size(), &energies[0], &fwhms[0], &energy_errors[0], &fwhm_errors[0]); }
}

// Fill accuracy histograms for reconstructed Eprot energy and theta
void FillAccuracyHistograms(TFile* file, const TString& prefix, const TString& suffix, TH2D* accuracyHist, const TString& sumSuffix = "", TH2D* summedHist = nullptr) {
    std::map<std::string, TKey*> latestKeys;

    // Find latest cycle for each object
    TIter next(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        std::string baseName = key->GetName();
        if (latestKeys.find(baseName) == latestKeys.end() || key->GetCycle() > latestKeys[baseName]->GetCycle())
            latestKeys[baseName] = key;
    }

    // Loop through latest keys
    for (const auto& pair : latestKeys) {
        TKey* key = pair.second;        TObject* obj = key->ReadObj();
        // --- Case 1: TH1D stitching vs excitation energy ---
        if (obj->InheritsFrom("TH1D")) {
            TH1D* hist = (TH1D*)obj;        TString name = hist->GetName();
            if (!name.BeginsWith(prefix)) continue;     if (!name.EndsWith(suffix)) continue;

            // Extract excitation energy from histogram name
            TString nameCopy = name;
            nameCopy.ReplaceAll(prefix, "");            nameCopy.ReplaceAll(suffix, "");            nameCopy.ReplaceAll("p", ".");
            double Eexc = nameCopy.Atof();

            for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
                double yValue  = hist->GetBinCenter(bin), content = hist->GetBinContent(bin), error = hist->GetBinError(bin);
                int xBin = accuracyHist->GetXaxis()->FindBin(Eexc), yBin = accuracyHist->GetYaxis()->FindBin(yValue);

                if (xBin >= 1 && xBin <= accuracyHist->GetNbinsX() && yBin >= 1 && yBin <= accuracyHist->GetNbinsY()) {
                    accuracyHist->SetBinContent(xBin, yBin, content);       accuracyHist->SetBinError(xBin, yBin, error);
                }
            }
        }
        // --- Case 2: TH2D summing (ignore excitation energy) ---
        else if (summedHist && obj->InheritsFrom("TH2D")) {
            TH2D* hist2D = (TH2D*)obj;      TString name = hist2D->GetName();
            if (!name.BeginsWith(prefix)) continue;      if (!name.EndsWith(sumSuffix)) continue;
            summedHist->Add(hist2D);
        }
    }
}


// Heavy residue coincidence plots in the XY plane
TH2D* histChannelPlots(const TString& reaction, const TString& prefix, const TString& suffix, TCanvas* canvas) {
    TH2D* sumHist = nullptr;
    std::map<double, int> colorMap;

    // Define channel labels
    std::vector<const char*> hrLabels = {"HRg", "HR1n", "HR2n", "HR3n"};

    // Color gradient definitions for each channel
    std::vector<std::vector<float>> gradients = {
        {0.0, 0.0, 0.5, 0.0, 1.0, 1.0}, // Blue to cyan
        {0.4, 0.0, 0.0, 1.0, 0.8, 0.0}, // Dark red to orange
        {0.0, 0.4, 0.0, 0.6, 1.0, 0.0}, // Green to lime
        {0.4, 0.0, 0.5, 1.0, 0.4, 0.7}  // Purple to pink
    };

    canvas->cd();
    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);

    for (size_t chanIdx = 0; chanIdx < hrLabels.size(); ++chanIdx) {
        const char* channel = hrLabels[chanIdx];
        TString fname = Form("../%s_sim/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction.Data(), reaction.Data(), channel);
        TFile* file = TFile::Open(fname);
        if (!file || file->IsZombie()) continue;

        // Collect latest keys
        std::map<std::string, TKey*> latestKeys;
        TIter next(file->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            std::string baseName = key->GetName();
            if (latestKeys.find(baseName) == latestKeys.end() || key->GetCycle() > latestKeys[baseName]->GetCycle()) {
                latestKeys[baseName] = key;
            }
        }

        // Collect histograms by excitation energy
        std::map<double, TH2D*> histogramsByEexc;
        for (auto& it : latestKeys) {
            TKey* latestKey = it.second;
            TObject* obj = latestKey->ReadObj();
            if (!obj || !obj->InheritsFrom("TH2D")) continue;
            TH2D* hist = dynamic_cast<TH2D*>(obj);
            if (!hist) continue;
            hist->SetDirectory(0);
            TString name = hist->GetName();
            if (!name.BeginsWith(prefix) || !name.EndsWith(suffix)) continue;
            TString nameCopy = name;
            nameCopy.ReplaceAll(prefix, "");    nameCopy.ReplaceAll(suffix, "");    nameCopy.ReplaceAll("p", ".");
            double Eexc = nameCopy.Atof();
            histogramsByEexc[Eexc] = hist;
        }

        // Color map for this channel
        std::vector<double> excEnergies;
        for (const auto& pair : histogramsByEexc) excEnergies.push_back(pair.first);
        const std::vector<float>& grad = gradients[chanIdx];
        float r1 = grad[0], g1 = grad[1], b1 = grad[2];
        float r2 = grad[3], g2 = grad[4], b2 = grad[5];
        colorMap = MapEnergiesToColors(excEnergies, r1, g1, b1, r2, g2, b2);

        // Plot and sum
        for (const auto& pair : histogramsByEexc) {
            double Eexc = pair.first;
            TH2D* hist = pair.second;
            
            if (!sumHist) {
                sumHist = (TH2D*)hist->Clone(Form("%s%s_sum", prefix.Data(), suffix.Data()));
                sumHist->SetDirectory(0);
            } else {
                sumHist->Add(hist);
            }
            hist->SetMarkerStyle(20);       hist->SetMarkerSize(0.3);       hist->SetMarkerColor(colorMap[Eexc]);
            
            if (canvas->GetListOfPrimitives()->GetSize() == 0) { hist->Draw("SCAT"); } 
            else { hist->Draw("SCAT SAME"); }
        }
        file->Close();
    }

    if (sumHist) {
        sumHist->SetTitle(prefix + suffix + " Summed Histogram over E*");
        canvas->Update();
        return sumHist;
    } else {
        std::cerr << "No matching histograms found for prefix " << prefix << std::endl;
        return nullptr;
    }
}

// Heavy residue coincidence survival fraction plots
void hrFractionPlot(TCanvas* canvas, const TString& prefix) {
    std::vector<TString> reactions = {"206Pbpp", "206Pbdd", "206Pbdp"}; 
    std::vector<const char*> hrLabels = {"HRg", "HR1n", "HR2n", "HR3n"};
    std::vector<int> markerStyles = {20, 21, 22};

    std::vector<std::vector<float>> gradients = {
        {0.0, 0.0, 0.5, 0.0, 1.0, 1.0},  // Blue to Cyan
        {0.4, 0.0, 0.0, 1.0, 0.8, 0.0},  // Dark Red to Orange
        {0.0, 0.4, 0.0, 0.6, 1.0, 0.0},  // Green to Lime
        {0.4, 0.0, 0.5, 1.0, 0.4, 0.7}   // Purple to Pink
    };

    canvas->cd();
    TMultiGraph* mg = new TMultiGraph();
    TLegend* legend = new TLegend(0.7, 0.47, 0.9, 0.67);    legend->SetBorderSize(0);       legend->SetFillStyle(0);
    TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);     legend2->SetBorderSize(0);      legend2->SetFillStyle(0);

    for (size_t chanIdx = 0; chanIdx < hrLabels.size(); ++chanIdx) {
        const char* channel = hrLabels[chanIdx];
        bool addedToLegend = false;  // one per channel
        std::vector<float> grad = gradients[chanIdx];
        float r1 = grad[0], g1 = grad[1], b1 = grad[2];
        float r2 = grad[3], g2 = grad[4], b2 = grad[5];

        for (size_t rx = 0; rx < reactions.size(); ++rx) {
            const TString& reaction = reactions[rx];
            TString fname = Form("../%s_sim/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction.Data(), reaction.Data(), channel);
            TFile* file = TFile::Open(fname);
            if (!file || file->IsZombie()) continue;

            std::map<std::string, TKey*> latestKeys;
            TIter next(file->GetListOfKeys());
            while (TKey* key = (TKey*)next()) {
                std::string baseName = key->GetName();
                if (latestKeys.find(baseName) == latestKeys.end() || key->GetCycle() > latestKeys[baseName]->GetCycle())
                    latestKeys[baseName] = key;
            }

            std::vector<std::pair<double, double>> points;

            for (std::map<std::string, TKey*>::iterator it = latestKeys.begin(); it != latestKeys.end(); ++it) {
                std::string nameStr = it->first;
                if (!TString(nameStr).BeginsWith(prefix) || !TString(nameStr).EndsWith("_HR_coinc")) continue;

                TKey* keyHR = it->second;
                TH2D* hrHist = dynamic_cast<TH2D*>(keyHR->ReadObj());
                if (!hrHist) continue;

                TString dipoleName = TString(nameStr).ReplaceAll("_HR_coinc", "_Dipole_coinc");
                TH2D* dipoleHist = dynamic_cast<TH2D*>(file->Get(dipoleName));
                if (!dipoleHist) continue;

                double hrEntries = hrHist->GetEntries();
                double dipoleEntries = dipoleHist->GetEntries();
                if (dipoleEntries <= 0) continue;

                TString energyStr = TString(nameStr.c_str());
                energyStr.ReplaceAll(prefix, "");
                energyStr.ReplaceAll("_HR_coinc", "");
                energyStr.ReplaceAll("p", ".");
                double Eexc = energyStr.Atof();
                double fraction = 100.0 * hrEntries / dipoleEntries;

                points.push_back(std::make_pair(Eexc, fraction));
            }

            std::sort(points.begin(), points.end());
            int n = points.size();
            if (n == 0) continue;

            std::vector<double> x(n), y(n);
            for (int i = 0; i < n; ++i) {
                x[i] = points[i].first;
                y[i] = points[i].second;
            }

            // Assign color per point based on excitation energy
            std::map<double, int> colorMap = MapEnergiesToColors(x, r1, g1, b1, r2, g2, b2);
            for (int i = 0; i < n; ++i) {
                TGraph* ptGraph = new TGraph(1, &x[i], &y[i]);
                ptGraph->SetMarkerStyle(markerStyles[rx]);
                ptGraph->SetMarkerSize(1.2);
                ptGraph->SetMarkerColor(colorMap[x[i]]);
                ptGraph->SetLineColor(colorMap[x[i]]);
                mg->Add(ptGraph);
            }

            if (!addedToLegend) {
                double midE = x[n / 2];     int midColor = colorMap[midE];
                TGraph* dummy = new TGraph();       dummy->SetMarkerStyle(20);      dummy->SetMarkerSize(1.2);
                dummy->SetMarkerColor(midColor);    dummy->SetLineColor(midColor);  dummy->SetPoint(0, 0, 0);
                legend->AddEntry(dummy, channel, "p");
                addedToLegend = true;
            }

        }
    }

    mg->Draw("AP");     mg->SetTitle(Form("Fraction of HR transmitted to detector (%s)", prefix.Data()));
    mg->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    mg->GetYaxis()->SetTitle("HR coinc / Dipole coinc (%)");

    std::map<TString, int> reactionMarker;
    reactionMarker["206Pbpp"] = 20;     reactionMarker["206Pbdd"] = 21;     reactionMarker["206Pbdp"] = 22;

    std::map<TString, TString> reactionLabel;
    reactionLabel["206Pbpp"] = "206Pb(p,p')";   reactionLabel["206Pbdd"] = "206Pb(d,d')";  reactionLabel["206Pbdp"] = "206Pb(d,p)";

    for (std::map<TString, int>::iterator it = reactionMarker.begin(); it != reactionMarker.end(); ++it) {
        const TString& reaction = it->first;    int markerStyle = it->second;
        TGraph* dummy = new TGraph();       dummy->SetMarkerStyle(markerStyle);     dummy->SetMarkerSize(1.2);       
        dummy->SetMarkerColor(kBlack);      dummy->SetLineColor(kBlack);            dummy->SetPoint(0, 0, 0); // dummy point
        legend2->AddEntry(dummy, reactionLabel[reaction], "p");
    }

    legend->Draw();
    legend2->Draw();
    canvas->Update();
}

// escaping particles from the E1 and Eres detectors: veto plots for dE position map
void vetoPlots(TCanvas* canvas, const TString& prefix, const TString& reaction) {
    std::vector<const char*> hrLabels = {"HRg", "HR1n", "HR2n", "HR3n"};
    
    canvas->cd();   canvas->Clear();
    
    // Create detector box (122 x 40 mm centered at 0,0)
    TBox* detectorBox = new TBox(-61, -20, 61, 20);  // 122x40 mm centered at (0,0)
    detectorBox->SetLineColor(kBlack);      detectorBox->SetLineWidth(2);       detectorBox->SetFillStyle(0);  // Transparent fill
    
    // Color definitions for E1 and Eres veto events
    int e1Color = kRed-3, eresColor = kAzure+2;
    
    // Collect all E1 and Eres veto histograms across all channels
    std::vector<TH2D*> e1VetoHists, eresVetoHists;
    
    for (size_t chanIdx = 0; chanIdx < hrLabels.size(); ++chanIdx) {
        const char* channel = hrLabels[chanIdx];
        TString fname = Form("../%s_sim/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction.Data(), reaction.Data(), channel);
        TFile* file = TFile::Open(fname);
        if (!file || file->IsZombie()) continue;

        // Collect latest keys
        std::map<std::string, TKey*> latestKeys;
        TIter next(file->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            std::string baseName = key->GetName();
            if (latestKeys.find(baseName) == latestKeys.end() || key->GetCycle() > latestKeys[baseName]->GetCycle()) {
                latestKeys[baseName] = key;
            }
        }
        
        // Find E1 and Eres veto histograms
        for (auto& it : latestKeys) {
            TKey* latestKey = it.second;
            TObject* obj = latestKey->ReadObj();
            if (!obj || !obj->InheritsFrom("TH2D")) continue;
            TH2D* hist = dynamic_cast<TH2D*>(obj);
            if (!hist) continue;
            hist->SetDirectory(0);
            TString name = hist->GetName();
            
            // Check if it's a veto histogram for this prefix
            if (name.BeginsWith(prefix) && name.EndsWith("_E1Veto_dE_posMap")) {
                e1VetoHists.push_back(hist);
            } else if (name.BeginsWith(prefix) && name.EndsWith("_EresVeto_dE_posMap")) {
                eresVetoHists.push_back(hist);
            }
        }
        file->Close();
    }
    
    // Create summed histograms
    TH2D* summedE1Veto = nullptr;    TH2D* summedEresVeto = nullptr;
    
    // Sum E1 veto histograms
    for (auto* hist : e1VetoHists) {
        if (!summedE1Veto) {
            summedE1Veto = (TH2D*)hist->Clone(Form("%s_E1Veto_sum", prefix.Data()));
            summedE1Veto->SetDirectory(0);
        } else {
            summedE1Veto->Add(hist);
        }
    }
    
    // Sum Eres veto histograms
    for (auto* hist : eresVetoHists) {
        if (!summedEresVeto) {
            summedEresVeto = (TH2D*)hist->Clone(Form("%s_EresVeto_sum", prefix.Data()));
            summedEresVeto->SetDirectory(0);
        } else {
            summedEresVeto->Add(hist);
        }
    }
    
    // Draw the histograms as scatter plots
    if (summedE1Veto) {
        summedE1Veto->SetMarkerStyle(20);       summedE1Veto->SetMarkerSize(0.3);       summedE1Veto->SetMarkerColor(e1Color);
        summedE1Veto->Draw("SCAT");
    }
    
    if (summedEresVeto) {
        summedEresVeto->SetMarkerStyle(21);     summedEresVeto->SetMarkerSize(0.3);     summedEresVeto->SetMarkerColor(eresColor);
        if (summedE1Veto) {
            summedEresVeto->Draw("SCAT SAME");
        } else {
            summedEresVeto->Draw("SCAT");
        }
    }
    
    // Draw detector box
    detectorBox->Draw("SAME");
    
    // Create legend - horizontal across the top outside the detector box
    TLegend* legend = new TLegend(0.15, 0.85, 0.9, 0.88);
    legend->SetBorderSize(0);   legend->SetFillStyle(0);    legend->SetTextSize(0.02);
    legend->SetNColumns(3);  // Arrange entries horizontally
    
    if (summedE1Veto) {
        TMarker* dummyE1 = new TMarker(0, 0, 20);
        dummyE1->SetMarkerColor(e1Color);       dummyE1->SetMarkerSize(1.0);
        legend->AddEntry(dummyE1, Form("E1 Veto Events (%d)", (int)summedE1Veto->GetEntries()), "p");
    }
    
    if (summedEresVeto) {
        TMarker* dummyEres = new TMarker(0, 0, 21);
        dummyEres->SetMarkerColor(eresColor);   dummyEres->SetMarkerSize(1.0);
        legend->AddEntry(dummyEres, Form("Eres Veto Events (%d)", (int)summedEresVeto->GetEntries()), "p");
    }
    
    legend->AddEntry(detectorBox, "Detector (122x40 mm)", "l");
    legend->Draw();
    
    // Set axis ranges to show the detector area clearly
    if (summedE1Veto || summedEresVeto) {
        TH2D* baseHist = summedE1Veto ? summedE1Veto : summedEresVeto;
        baseHist->GetXaxis()->SetRangeUser(-80, 80);
        baseHist->GetYaxis()->SetRangeUser(-30, 30);
        baseHist->SetTitle(Form("%s Veto Position Map; X Position (mm); Y Position (mm)", prefix.Data()));
    }
    
    // Disable gStats display
    gStyle->SetOptStat(0);
    
    canvas->Update();
}

// Function to create veto-gated banana plots overlaid on ungated bananas
void vetoBananaPlots(TCanvas* canvasE1, TCanvas* canvasEres, const TString& prefix, const TString& reaction) {
    std::vector<const char*> hrLabels = {"HRg", "HR1n", "HR2n", "HR3n"};
    
    // Collect all ungated and veto-gated banana histograms across all channels
    std::vector<TH2D*> ungatedE1Hists, ungatedEresHists, vetoE1Hists, vetoEresHists;
    
    for (size_t chanIdx = 0; chanIdx < hrLabels.size(); ++chanIdx) {
        const char* channel = hrLabels[chanIdx];
        TString fname = Form("../%s_sim/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction.Data(), reaction.Data(), channel);
        TFile* file = TFile::Open(fname);
        if (!file || file->IsZombie()) continue;

        // Collect latest keys
        std::map<std::string, TKey*> latestKeys;
        TIter next(file->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            std::string baseName = key->GetName();
            if (latestKeys.find(baseName) == latestKeys.end() || key->GetCycle() > latestKeys[baseName]->GetCycle()) {
                latestKeys[baseName] = key;
            }
        }
        
        // Find ungated and veto-gated banana histograms
        for (auto& it : latestKeys) {
            TKey* latestKey = it.second;
            TObject* obj = latestKey->ReadObj();
            if (!obj || !obj->InheritsFrom("TH2D")) continue;
            TH2D* hist = dynamic_cast<TH2D*>(obj);
            if (!hist) continue;
            hist->SetDirectory(0);
            TString name = hist->GetName();
            
            // Check for different histogram types
            if (name.BeginsWith(prefix) && name.EndsWith("_dE_v_E1_Banana") && !name.Contains("_veto")) {
                ungatedE1Hists.push_back(hist);
            } else if (name.BeginsWith(prefix) && name.EndsWith("_dE_v_Eres_Banana") && !name.Contains("_veto")) {
                ungatedEresHists.push_back(hist);
            } else if (name.BeginsWith(prefix) && name.EndsWith("_dE_v_E1_Banana_veto")) {
                vetoE1Hists.push_back(hist);
            } else if (name.BeginsWith(prefix) && name.EndsWith("_dE_v_Eres_Banana_veto")) {
                vetoEresHists.push_back(hist);
            }
        }
        file->Close();
    }
    
    // Create summed histograms
    TH2D* summedUngatedE1 = nullptr, *summedUngatedEres = nullptr, *summedVetoE1 = nullptr, *summedVetoEres = nullptr;
    
    for (auto* hist : ungatedE1Hists) {     // Sum ungated E1 histograms
        if (!summedUngatedE1) { summedUngatedE1 = (TH2D*)hist->Clone(Form("%s_ungated_E1_sum", prefix.Data()));         summedUngatedE1->SetDirectory(0);   }
        else { summedUngatedE1->Add(hist); }
    }
    for (auto* hist : ungatedEresHists) {   // Sum ungated Eres histograms
        if (!summedUngatedEres) { summedUngatedEres = (TH2D*)hist->Clone(Form("%s_ungated_Eres_sum", prefix.Data()));         summedUngatedEres->SetDirectory(0);   }
        else { summedUngatedEres->Add(hist); }
    }
    for (auto* hist : vetoE1Hists) {        // Sum veto E1 histograms
        if (!summedVetoE1) { summedVetoE1 = (TH2D*)hist->Clone(Form("%s_veto_E1_sum", prefix.Data()));         summedVetoE1->SetDirectory(0);   }
        else { summedVetoE1->Add(hist); }
    }
    for (auto* hist : vetoEresHists) {      // Sum veto Eres histograms
        if (!summedVetoEres) { summedVetoEres = (TH2D*)hist->Clone(Form("%s_veto_Eres_sum", prefix.Data()));         summedVetoEres->SetDirectory(0);   }
        else { summedVetoEres->Add(hist); }
    }
    
    gStyle->SetPalette(kDeepSea);       // Set palette for ungated histograms

    // Draw E1 banana plot
    canvasE1->cd();     canvasE1->Clear();
    if (summedUngatedE1) { 
        summedUngatedE1->Draw("COLZ");     summedUngatedE1->SetMaximum(50);
        summedUngatedE1->SetTitle(Form("%s Veto-Gated E1 Banana Plot; dE (MeV); E1 (MeV)", prefix.Data()));
    }
    if (summedVetoE1) {
        summedVetoE1->SetMarkerStyle(20);       summedVetoE1->SetMarkerSize(0.3);       summedVetoE1->SetMarkerColor(kRed);
        summedVetoE1->Draw("SCAT SAME");
    }
    
    // Create legend for E1 plot
    TLegend* legendE1 = new TLegend(0.55, 0.8, 0.75, 0.9);
    legendE1->SetBorderSize(0);   legendE1->SetFillStyle(0);    legendE1->SetTextSize(0.03);
    
    if (summedUngatedE1) {
        TMarker* dummyUngatedE1 = new TMarker(0, 0, 20);
        dummyUngatedE1->SetMarkerColor(kCyan);    dummyUngatedE1->SetMarkerSize(1.0);
        legendE1->AddEntry(dummyUngatedE1, Form("Ungated E1 Events (%d)", (int)summedUngatedE1->GetEntries()), "p");
    }
    if (summedVetoE1) {
        TMarker* dummyVetoE1 = new TMarker(0, 0, 20);
        dummyVetoE1->SetMarkerColor(kRed);        dummyVetoE1->SetMarkerSize(1.0);
        legendE1->AddEntry(dummyVetoE1, Form("Veto-Gated E1 Events (%d)", (int)summedVetoE1->GetEntries()), "p");
    }
    legendE1->Draw();
    gStyle->SetOptStat(0);
    canvasE1->Update();
    
    // Draw Eres banana plot
    canvasEres->cd();
    canvasEres->Clear();
    if (summedUngatedEres) { 
        summedUngatedEres->Draw("COLZ");     summedUngatedEres->SetMaximum(50);
        summedUngatedEres->SetTitle(Form("%s Veto-Gated Eres Banana Plot; dE (MeV); Eres (MeV)", prefix.Data()));
    }
    if (summedVetoEres) {
        summedVetoEres->SetMarkerStyle(21);     summedVetoEres->SetMarkerSize(0.3);     summedVetoEres->SetMarkerColor(kRed);
        summedVetoEres->Draw("SCAT SAME");
    }
    
    // Create legend for Eres plot
    TLegend* legendEres = new TLegend(0.55, 0.8, 0.75, 0.9);
    legendEres->SetBorderSize(0);   legendEres->SetFillStyle(0);    legendEres->SetTextSize(0.03);
    
    if (summedUngatedEres) {
        TMarker* dummyUngatedEres = new TMarker(0, 0, 21);
        dummyUngatedEres->SetMarkerColor(kCyan);    dummyUngatedEres->SetMarkerSize(1.0);
        legendEres->AddEntry(dummyUngatedEres, Form("Ungated Eres Events (%d)", (int)summedUngatedEres->GetEntries()), "p");
    }
    if (summedVetoEres) {
        TMarker* dummyVetoEres = new TMarker(0, 0, 21);
        dummyVetoEres->SetMarkerColor(kRed);        dummyVetoEres->SetMarkerSize(1.0);
        legendEres->AddEntry(dummyVetoEres, Form("Veto-Gated Eres Events (%d)", (int)summedVetoEres->GetEntries()), "p");
    }
    legendEres->Draw();
    gStyle->SetOptStat(0);
    canvasEres->Update();
} 