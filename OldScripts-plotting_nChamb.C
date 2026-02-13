// --
// NECTAR Excitation energy resolution: fits excitation energy histograms and plots
// Last modification : G. Leckenby - 1/07/2025

// ---------------------------------------------------------------------------------------
// ################## Define which plots you want to write ###############################
// ---------------------------------------------------------------------------------------
const bool write_excRes = false;            // Excitation energy resolution graph (define filelist for additional series)
const bool write_accuracy = true;          // Reconstructed Eprot energy and theta accuracy histograms
const bool write_full_banana = false;       // Full banana plots (dE v E1 and dE v Eres) summed over all angles
const bool write_banana_vs_theta = false;    // Banana plots as above as a function of theta (need to define theta histograms)
const bool write_HR_XY_plots = false;       // Heavy residue coincidence plots in the XY plane
const bool write_HR_X_Eexc = false;          // Heavy residue coincidence plots with X position as a function of E*
const bool write_HR_frac = false;           // Heavy residue transmission graphs for the 3 different reactions.
const bool write_veto_plots = false;         // Veto position maps with overlayed E1 and Eres events
// ---------------------------------------------------------------------------------------

// Include the plotting functions header and implementation
#include "UtilityScripts/funcs_plotting.h"
#include "UtilityScripts/funcs_plotting.C"


// main execution function
void plotting_nChamb() {
    // Excitation energy resolution plot
    if(write_excRes) {
        std::vector<TString> reactions = {"206Pbdp", "206Pbdd"};
        std::vector<TString> hrLabels = {"HRg", "HR1n", "HR2n", "HR3n"};
        std::vector<TGraphErrors*> vecGraphs;
        TCanvas* c1a = new TCanvas("c1a", "E* Res (sigma) vs Eexc", 50, 50, 800, 600);
        
        // Create separate resolution histograms for each reaction + FP/BP combo
        TH2D* resolutionHistFP_dp = new TH2D("resolution_histFP_dp", "FP 206Pb(d,p) Reconstructed E* Values; E* (MeV); Ereco - E* (MeV)", 26, 0.0, 26.0, 100, -5.0, 5.0);
        TH2D* resolutionHistBP_dp = new TH2D("resolution_histBP_dp", "BP 206Pb(d,p) Reconstructed E* Values; E* (MeV); Ereco - E* (MeV)", 26, 0.0, 26.0, 100, -5.0, 5.0);
        TH2D* resolutionHistFP_dd = new TH2D("resolution_histFP_dd", "FP 206Pb(d,d') Reconstructed E* Values; E* (MeV); Ereco - E* (MeV)", 26, 0.0, 26.0, 100, -5.0, 5.0);
        TH2D* resolutionHistBP_dd = new TH2D("resolution_histBP_dd", "BP 206Pb(d,d') Reconstructed E* Values; E* (MeV); Ereco - E* (MeV)", 26, 0.0, 26.0, 100, -5.0, 5.0);

        for (size_t reactIdx = 0; reactIdx < reactions.size(); ++reactIdx) {        // Process each reaction
            TString reaction = reactions[reactIdx];
            
            for (size_t chanIdx = 0; chanIdx < hrLabels.size(); ++chanIdx) {        // Process each channel for this reaction
                TString channel = hrLabels[chanIdx];
                TString reaction_dir = getReactionDir(reaction);
                TString fname = Form("../%s/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction_dir.Data(), reaction.Data(), channel.Data());
                TFile* file = TFile::Open(fname);
                if (!file || file->IsZombie()) continue;

                TH2D* histFP = (reactIdx == 0) ? resolutionHistFP_dp : resolutionHistFP_dd;     // Choose appropriate resolution histograms based on reaction
                TH2D* histBP = (reactIdx == 0) ? resolutionHistBP_dp : resolutionHistBP_dd;
                TGraphErrors* graphFP = FitResolutionGraph(file, "h_FP_Eex", histFP);       if (graphFP) vecGraphs.push_back(graphFP);
                TGraphErrors* graphBP = FitResolutionGraph(file, "h_BP_Eex", histBP);       if (graphBP) vecGraphs.push_back(graphBP);
                file->Close();
            }
        }

        std::vector<int> colors = {kAzure+7, kAzure+4, kCyan-3, kTeal-1, kRed+1, kPink+8, kMagenta+2, kMagenta-9};
        std::vector<int> markers = {20, 24};
        TLegend* legend = new TLegend(0.7, 0.6, 0.9, 0.9);
        TMultiGraph* mg = new TMultiGraph("mg", "E* Resolution (sigma)");
        
        for (size_t i = 0; i < vecGraphs.size(); ++i) {
            auto* graph = vecGraphs[i];
            graph->SetMarkerStyle(markers[i % 2]);      graph->SetMarkerColor(colors[i/2]);
            graph->SetMarkerSize(0.7);                  graph->SetLineColor(colors[i/2]);
            mg->Add(graph, "P");
            
            int reactIdx = i / 8;       int chanIdx = (i % 8) / 2;       int telIdx = i % 2;  
            TString legendEntry = Form("%s %s %s", reactions[reactIdx].Data(), hrLabels[chanIdx].Data(), (telIdx == 0) ? "FP" : "BP");
            legend->AddEntry(graph, legendEntry, "p");
        }
        mg->Draw("A");      mg->GetYaxis()->SetRangeUser(0.0, 1.0);
        legend->Draw();     c1a->Update();
        
        // Plot separate resolution histograms for each reaction + telescope combo
        TCanvas* c1b = new TCanvas("c1b", "FP 206Pb(d,p) Resolution Map", 900, 50, 800, 600);
        resolutionHistFP_dp->Draw("COLZ");
        TLine* zeroLine = new TLine(0.0, 0.0, 26.0, 0.0);     zeroLine->SetLineColor(kRed);     zeroLine->SetLineWidth(2);
        zeroLine->Draw("SAME");        
        TCanvas* c1c = new TCanvas("c1c", "BP 206Pb(d,p) Resolution Map", 900, 700, 800, 600);
        resolutionHistBP_dp->Draw("COLZ");      zeroLine->Draw("SAME");        
        TCanvas* c1d = new TCanvas("c1d", "FP 206Pb(d,d') Resolution Map", 1750, 50, 800, 600);
        resolutionHistFP_dd->Draw("COLZ");      zeroLine->Draw("SAME");
        TCanvas* c1e = new TCanvas("c1e", "BP 206Pb(d,d') Resolution Map", 1750, 700, 800, 600);
        resolutionHistBP_dd->Draw("COLZ");      zeroLine->Draw("SAME");
    }
    
    // Reconstructed Eprot energy and theta accuracy plots
    if(write_accuracy) {
        TString reaction = "206Pbdp";
        std::vector<TString> hrLabels = {"HRg", "HR1n", "HR2n", "HR3n"};
        
        // Create separate accuracy histograms for each reaction + FP/BP combo
        TH2D* accuracyHistFP_dp_Eprot = new TH2D("accuracy_histFP_dp_Eprot", "FP 206Pb(d,p) Reconstructed Eprot Values; E* (MeV); Ereco_prot - Etrue_prot (MeV)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dp_Eprot = new TH2D("accuracy_histBP_dp_Eprot", "BP 206Pb(d,p) Reconstructed Eprot Values; E* (MeV); Ereco_prot - Etrue_prot (MeV)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        TH2D* accuracyHistFP_dd_Eprot = new TH2D("accuracy_histFP_dd_Eprot", "FP 206Pb(d,d') Reconstructed Eprot Values; E* (MeV); Ereco_prot - Etrue_prot (MeV)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dd_Eprot = new TH2D("accuracy_histBP_dd_Eprot", "BP 206Pb(d,d') Reconstructed Eprot Values; E* (MeV); Ereco_prot - Etrue_prot (MeV)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        
        TH2D* accuracyHistFP_dp_EprotvEprot = new TH2D("accuracy_histFP_dp_EprotvEprot", "FP 206Pb(d,p) Reconstructed Eprot Values vs Eprot; Eprot (MeV); Ereco_prot - Etrue_prot (MeV)", 300, 0, 150, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dp_EprotvEprot = new TH2D("accuracy_histBP_dp_EprotvEprot", "BP 206Pb(d,p) Reconstructed Eprot Values vs Eprot; Eprot (MeV); Ereco_prot - Etrue_prot (MeV)", 300, 0, 150, 200, -5.0, 5.0);
        TH2D* accuracyHistFP_dd_EprotvEprot = new TH2D("accuracy_histFP_dd_EprotvEprot", "FP 206Pb(d,d') Reconstructed Eprot Values vs Eprot; Eprot (MeV); Ereco_prot - Etrue_prot (MeV)", 300, 0, 150, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dd_EprotvEprot = new TH2D("accuracy_histBP_dd_EprotvEprot", "BP 206Pb(d,d') Reconstructed Eprot Values vs Eprot; Eprot (MeV); Ereco_prot - Etrue_prot (MeV)", 300, 0, 150, 200, -5.0, 5.0);
        
        TH2D* accuracyHistFP_dp_theta = new TH2D("accuracy_histFP_dp_theta", "FP 206Pb(d,p) Reconstructed Theta Values; E* (MeV); Theta_reco - Theta_true (deg)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dp_theta = new TH2D("accuracy_histBP_dp_theta", "BP 206Pb(d,p) Reconstructed Theta Values; E* (MeV); Theta_reco - Theta_true (deg)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        TH2D* accuracyHistFP_dd_theta = new TH2D("accuracy_histFP_dd_theta", "FP 206Pb(d,d') Reconstructed Theta Values; E* (MeV); Theta_reco - Theta_true (deg)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dd_theta = new TH2D("accuracy_histBP_dd_theta", "BP 206Pb(d,d') Reconstructed Theta Values; E* (MeV); Theta_reco - Theta_true (deg)", 26, 0.0, 26.0, 200, -5.0, 5.0);
        
        TH2D* accuracyHistFP_dp_ThetavTheta = new TH2D("accuracy_histFP_dp_ThetavTheta", "FP 206Pb(d,p) Reconstructed Theta Values vs Theta; Theta (deg); Theta_reco - Theta_true (dge)", 180, 0, 90, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dp_ThetavTheta = new TH2D("accuracy_histBP_dp_ThetavTheta", "BP 206Pb(d,p) Reconstructed Theta Values vs Theta; Theta (deg); Theta_reco - Theta_true (dge)", 180, 0, 90, 200, -5.0, 5.0);
        TH2D* accuracyHistFP_dd_ThetavTheta = new TH2D("accuracy_histFP_dd_ThetavTheta", "FP 206Pb(d,d') Reconstructed Theta Values vs Theta; Theta (deg); Theta_reco - Theta_true (dge)", 180, 0, 90, 200, -5.0, 5.0);
        TH2D* accuracyHistBP_dd_ThetavTheta = new TH2D("accuracy_histBP_dd_ThetavTheta", "BP 206Pb(d,d') Reconstructed Theta Values vs Theta; Theta (deg); Theta_reco - Theta_true (dge)", 180, 0, 90, 200, -5.0, 5.0);

        for (size_t chanIdx = 0; chanIdx < hrLabels.size(); ++chanIdx) {
            TString channel = hrLabels[chanIdx];
            TString fname = Form("../%s_sim/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction.Data(), reaction.Data(), channel.Data());
            TFile* file = TFile::Open(fname);
            if (!file || file->IsZombie()) continue;
            // Choose appropriate accuracy histograms based on reaction
            TH2D* histFP_Eprot = (reaction == "206Pbdp") ? accuracyHistFP_dp_Eprot : accuracyHistFP_dd_Eprot;
            TH2D* histBP_Eprot = (reaction == "206Pbdp") ? accuracyHistBP_dp_Eprot : accuracyHistBP_dd_Eprot;
            TH2D* histFP_EprotvEprot = (reaction == "206Pbdp") ? accuracyHistFP_dp_EprotvEprot : accuracyHistFP_dd_EprotvEprot;
            TH2D* histBP_EprotvEprot = (reaction == "206Pbdp") ? accuracyHistBP_dp_EprotvEprot : accuracyHistBP_dd_EprotvEprot;
            TH2D* histFP_theta = (reaction == "206Pbdp") ? accuracyHistFP_dp_theta : accuracyHistFP_dd_theta;
            TH2D* histBP_theta = (reaction == "206Pbdp") ? accuracyHistBP_dp_theta : accuracyHistBP_dd_theta;
            TH2D* histFP_ThetavTheta = (reaction == "206Pbdp") ? accuracyHistFP_dp_ThetavTheta : accuracyHistFP_dd_ThetavTheta;
            TH2D* histBP_ThetavTheta = (reaction == "206Pbdp") ? accuracyHistBP_dp_ThetavTheta : accuracyHistBP_dd_ThetavTheta;
            // Fill Eprot accuracy histograms
            FillAccuracyHistograms(file, "h_FP_Eex", "_Eprot_Accur", histFP_Eprot, "_Eprot_Accur_vs_Eprot", histFP_EprotvEprot);
            FillAccuracyHistograms(file, "h_BP_Eex", "_Eprot_Accur", histBP_Eprot, "_Eprot_Accur_vs_Eprot", histBP_EprotvEprot);
            // Fill theta accuracy histograms
            FillAccuracyHistograms(file, "h_FP_Eex", "_Theta_Accur", histFP_theta, "_Theta_Accur_vs_Theta", histFP_ThetavTheta);
            FillAccuracyHistograms(file, "h_BP_Eex", "_Theta_Accur", histBP_theta, "_Theta_Accur_vs_Theta", histBP_ThetavTheta);
            
            file->Close();
        }
        
        TLine* zeroLine = new TLine(0.0, 0.0, 26.0, 0.0);   zeroLine->SetLineColor(kRed);   zeroLine->SetLineWidth(2);
        TLine* zeroLineEprot = new TLine(0.0, 0.0, 150, 0.0);   zeroLineEprot->SetLineColor(kRed);   zeroLineEprot->SetLineWidth(2);
        TLine* zeroLineTheta = new TLine(0.0, 0.0, 90, 0.0);   zeroLineTheta->SetLineColor(kRed);   zeroLineTheta->SetLineWidth(2);           
        if (reaction == "206Pbdp") {
            // Plot (d,p) accuracy histograms
            TCanvas* c1f = new TCanvas("c1f", "FP 206Pb(d,p) Eprot Accuracy Map", 50, 50, 800, 600);
            accuracyHistFP_dp_Eprot->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1g = new TCanvas("c1g", "BP 206Pb(d,p) Eprot Accuracy Map", 50, 700, 800, 600);
            accuracyHistBP_dp_Eprot->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1h = new TCanvas("c1h", "FP 206Pb(d,p) Eprot Accuracy v Eprot", 150, 150, 800, 600);
            accuracyHistFP_dp_EprotvEprot->Draw("COLZ");    zeroLineEprot->Draw("SAME");
            TCanvas* c1i = new TCanvas("c1i", "BP 206Pb(d,p) Eprot Accuracy v Eprot", 150, 850, 800, 600);
            accuracyHistBP_dp_EprotvEprot->Draw("COLZ");    zeroLineEprot->Draw("SAME");
            TCanvas* c1j = new TCanvas("c1j", "FP 206Pb(d,p) Theta Accuracy Map", 900, 50, 800, 600);
            accuracyHistFP_dp_theta->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1k = new TCanvas("c1k", "BP 206Pb(d,p) Theta Accuracy Map", 900, 700, 800, 600);
            accuracyHistBP_dp_theta->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1l = new TCanvas("c1l", "FP 206Pb(d,p) Theta Accuracy v Theta", 1050, 150, 800, 600);
            accuracyHistFP_dp_ThetavTheta->Draw("COLZ");    zeroLineTheta->Draw("SAME");
            TCanvas* c1m = new TCanvas("c1m", "BP 206Pb(d,p) Theta Accuracy v Theta", 1050, 850, 800, 600);
            accuracyHistBP_dp_ThetavTheta->Draw("COLZ");    zeroLineTheta->Draw("SAME");
        } else {
            // Plot (d,d') accuracy histograms  
            TCanvas* c1f = new TCanvas("c1f", "FP 206Pb(d,d') Eprot Accuracy Map", 50, 50, 800, 600);
            accuracyHistFP_dd_Eprot->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1g = new TCanvas("c1g", "BP 206Pb(d,d') Eprot Accuracy Map", 50, 700, 800, 600);
            accuracyHistBP_dd_Eprot->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1h = new TCanvas("c1h", "FP 206Pb(d,d') Eprot Accuracy v Eprot", 150, 150, 800, 600);
            accuracyHistFP_dd_EprotvEprot->Draw("COLZ");    zeroLineEprot->Draw("SAME");
            TCanvas* c1i = new TCanvas("c1i", "BP 206Pb(d,d') Eprot Accuracy v Eprot", 150, 850, 800, 600);
            accuracyHistBP_dd_EprotvEprot->Draw("COLZ");    zeroLineEprot->Draw("SAME");
            TCanvas* c1j = new TCanvas("c1j", "FP 206Pb(d,d') Theta Accuracy Map", 900, 50, 800, 600);
            accuracyHistFP_dd_theta->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1k = new TCanvas("c1k", "BP 206Pb(d,d') Theta Accuracy Map", 900, 700, 800, 600);
            accuracyHistBP_dd_theta->Draw("COLZ");          zeroLine->Draw("SAME");
            TCanvas* c1l = new TCanvas("c1l", "FP 206Pb(d,d') Theta Accuracy v Theta", 1050, 150, 800, 600);
            accuracyHistFP_dd_ThetavTheta->Draw("COLZ");    zeroLineTheta->Draw("SAME");
            TCanvas* c1m = new TCanvas("c1m", "BP 206Pb(d,d') Theta Accuracy v Theta", 1050, 850, 800, 600);
            accuracyHistBP_dd_ThetavTheta->Draw("COLZ");    zeroLineTheta->Draw("SAME");
        }        
    }
    
    // Banana plots for telescopes 
    if(write_full_banana) {
        TString reaction1 = "206Pbdd";      TString reaction2 = "206Pbdp";
        TH2* h = nullptr;   TList* prims = nullptr;

        bool isFront = true;    TString prefix = "";    Float_t dEvE1XRangeHigh, dEvE1YRangeHigh, dEvEresXRangeHigh, dEvEresYRangeHigh;
        if (isFront) {  prefix = "h_FP_Eex";    dEvE1XRangeHigh = 25;   dEvE1YRangeHigh = 10;   dEvEresXRangeHigh = 70;     dEvEresYRangeHigh = 30; }
        else {          prefix = "h_BP_Eex";    dEvE1XRangeHigh = 12;   dEvE1YRangeHigh = 2.5;  dEvEresXRangeHigh = 100;    dEvEresYRangeHigh = 14; }
        
        TCanvas* c2a = new TCanvas("c2a", "dE v E1 Banana as a function of E*", 50, 50, 800, 600);
        TCanvas* c2b = new TCanvas("c2b", "dE v E1 Summed Banana", 900, 50, 800, 600);
        //TH2D* full_dEvE1_banana = nullptr;
        TH2D* full_dEvE1_banana_dd = histChannelPlots(reaction1, prefix, "_dE_v_E1_Banana", c2a);
        TH2D* full_dEvE1_banana_dp = histChannelPlots(reaction2, prefix, "_dE_v_E1_Banana", c2a);
        SetTH2RangeOnCanvas(c2a, 0, dEvE1XRangeHigh, 0, dEvE1YRangeHigh);

        TCanvas* c2c = new TCanvas("c2c", "dE+E1 v Eres Banana as a function of E*", 50, 700, 800, 600);
        TCanvas* c2d = new TCanvas("c2d", "dE+E1 v Eres Summed Banana", 900, 700, 800, 600);
        //TH2D* full_dEvE1_banana = nullptr;
        TH2D* full_dEvEres_banana_dd = histChannelPlots(reaction1, prefix, "_dE_v_Eres_Banana", c2c);
        TH2D* full_dEvEres_banana_dp = histChannelPlots(reaction2, prefix, "_dE_v_Eres_Banana", c2c);
        SetTH2RangeOnCanvas(c2c, 0, dEvEresXRangeHigh, 0, dEvEresYRangeHigh);

        Int_t maxVal = 50;
        if (isFront) { gStyle->SetPalette(kLake); } else { gStyle->SetPalette(kSunset); maxVal = 100; }
        full_dEvE1_banana_dd->SetMaximum(maxVal);      full_dEvE1_banana_dp->SetMaximum(maxVal); 
        full_dEvE1_banana_dd->GetXaxis()->SetRangeUser(0, dEvE1XRangeHigh);          
        full_dEvE1_banana_dd->GetYaxis()->SetRangeUser(0, dEvE1YRangeHigh);
        full_dEvEres_banana_dd->SetMaximum(maxVal);    full_dEvEres_banana_dp->SetMaximum(maxVal);
        full_dEvEres_banana_dd->GetXaxis()->SetRangeUser(0, dEvEresXRangeHigh);       
        full_dEvEres_banana_dd->GetYaxis()->SetRangeUser(0, dEvEresYRangeHigh);
        // Add legends with entry counts for dE v E1 banana
        c2b->cd();       full_dEvE1_banana_dd->Draw("COLZ");     full_dEvE1_banana_dp->Draw("COLZ SAME");
        
        TLegend* legend_dEvE1 = new TLegend(0.6, 0.8, 0.8, 0.9);
        legend_dEvE1->SetBorderSize(0);     legend_dEvE1->SetFillStyle(0);      legend_dEvE1->SetTextSize(0.03);
        
        // Create dummy histograms for legend entries
        TH2D* dummy_dd = new TH2D("dummy_dd", "", 1, 0, 1, 1, 0, 1);
        dummy_dd->SetFillColor(kRed);       dummy_dd->SetLineColor(kRed);
        
        TH2D* dummy_dp = new TH2D("dummy_dp", "", 1, 0, 1, 1, 0, 1);
        dummy_dp->SetFillColor(kBlue);      dummy_dp->SetLineColor(kBlue);
        
        legend_dEvE1->AddEntry(dummy_dd, Form("206Pb(d,d') (%d)", (int)full_dEvE1_banana_dd->GetEntries()), "f");
        legend_dEvE1->AddEntry(dummy_dp, Form("206Pb(d,p) (%d)", (int)full_dEvE1_banana_dp->GetEntries()), "f");
        legend_dEvE1->Draw();
        c2b->Update();
        
        // Add legends with entry counts for dE v Eres banana
        c2d->cd();       full_dEvEres_banana_dd->Draw("COLZ");   full_dEvEres_banana_dp->Draw("COLZ SAME");
        
        TLegend* legend_dEvEres = new TLegend(0.6, 0.8, 0.8, 0.9);
        legend_dEvEres->SetBorderSize(0);   legend_dEvEres->SetFillStyle(0);    legend_dEvEres->SetTextSize(0.03);
        
        legend_dEvEres->AddEntry(dummy_dd, Form("206Pb(d,d') (%d)", (int)full_dEvEres_banana_dd->GetEntries()), "f");
        legend_dEvEres->AddEntry(dummy_dp, Form("206Pb(d,p) (%d)", (int)full_dEvEres_banana_dp->GetEntries()), "f");
        legend_dEvEres->Draw();
        c2d->Update();
    }

    // Banana plots as a function of theta
    if(write_banana_vs_theta) {
        TString reaction1 = "206Pbdd";      TString reaction2 = "206Pbdp";
        
        TCanvas* c3a = new TCanvas("c3a", "dE v E1 Banana as a function of E* (70 < #theta < 80)", 900, 50, 800, 600);
        TCanvas* c3b = new TCanvas("c3b", "dE v E1 Summed Banana (70 < #theta < 80)", 50, 50, 800, 600);
        TH2D* theta1_dEvE1_banana_dd = histChannelPlots(reaction1, "h_Eex", "_dE_v_E1_Banana_theta_70_80", c3a);
        TH2D* theta1_dEvE1_banana_dp = histChannelPlots(reaction2, "h_Eex", "_dE_v_E1_Banana_theta_70_80", c3a);
        theta1_dEvE1_banana_dd->SetMaximum(20);     theta1_dEvE1_banana_dp->SetMaximum(20);
        c3a->Update();

        TCanvas* c3c = new TCanvas("c3c", "dE v Eres Banana as a function of E* (70 < #theta < 80)", 900, 700, 800, 600);
        TCanvas* c3d = new TCanvas("c3d", "dE v Eres Summed Banana (70 < #theta < 80)", 50, 700, 800, 600);
        TH2D* theta1_dEvEres_banana_dd = histChannelPlots(reaction1, "h_Eex", "_dE_v_ERes_Banana_theta_70_80", c3c);
        TH2D* theta1_dEvEres_banana_dp = histChannelPlots(reaction2, "h_Eex", "_dE_v_ERes_Banana_theta_70_80", c3c);
        theta1_dEvEres_banana_dd->SetMaximum(20);     theta1_dEvEres_banana_dp->SetMaximum(20);
        c3c->Update();

        gStyle->SetPalette(kLake);      
        c3b->cd();   theta1_dEvE1_banana_dd->Draw("COLZ");       theta1_dEvE1_banana_dp->Draw("COLZ SAME");
        c3d->cd();   theta1_dEvEres_banana_dd->Draw("COLZ");     theta1_dEvEres_banana_dp->Draw("COLZ SAME");
        c3b->Update();   c3d->Update();
    }
    
    // Heavy residue coincidence plots in the XY plane
    if(write_HR_XY_plots) {
        TBox* Dip_box = new TBox(-100, -35, 100, 35);
        Dip_box->SetLineColor(kRed);     Dip_box->SetLineWidth(3);    Dip_box->SetFillStyle(0);  // transparent
        TBox* HR_box = new TBox(-122-15, -20, -15, 20);
        HR_box->SetLineColor(kRed);     HR_box->SetLineWidth(3);    HR_box->SetFillStyle(0);  // transparent

        TString reaction = "206Pbdp";       
        TH2* h = nullptr;   TList* prims = nullptr;

        TCanvas* c4a = new TCanvas("c4a","FP Dipole entrance coincidence plot", 50, 50, 1000, 600);
        TH2D* h_FP_Sum_Dipole_coinc = histChannelPlots(reaction, "h_FP_Eex", "_Dipole_coinc", c4a);     
        Dip_box->Draw("SAME");
        SetTH2RangeOnCanvas(c4a, -150, 150, -50, 50);
        c4a->Update(); 
        
        TCanvas* c4b = new TCanvas("c4b","FP Heavy residue coincidence plot", 1100, 50, 1000, 600);
        TH2D* h_FP_Sum_HR_coinc = histChannelPlots(reaction, "h_FP_Eex", "_HR_coinc", c4b);        
        HR_box->Draw("SAME");
        SetTH2RangeOnCanvas(c4b, -150, 0, -30, 30);
        c4b->Update(); 

        TCanvas* c4c = new TCanvas("c4c","BP Dipole entrance coincidence plot", 50, 700, 1000, 600);
        TH2D* h_BP_Sum_Dipole_coinc = histChannelPlots(reaction, "h_BP_Eex", "_Dipole_coinc", c4c);     
        Dip_box->Draw("SAME");
        SetTH2RangeOnCanvas(c4c, -150, 150, -50, 50);
        c4c->Update(); 
        
        TCanvas* c4d = new TCanvas("c4d","BP Heavy residue coincidence plot", 1100, 700, 1000, 600);
        TH2D* h_BP_Sum_HR_coinc = histChannelPlots(reaction, "h_BP_Eex", "_HR_coinc", c4d);        
        HR_box->Draw("SAME");
        SetTH2RangeOnCanvas(c4d, -150, 0, -30, 30);
        c4d->Update();         
    }
    
    // Heavy residue coincidence plots with X position as a function of E*
    if(write_HR_X_Eexc) {
        TString reaction = "206Pbdp";
        TCanvas* c5a = new TCanvas("c5a","FP Heavy residue X vs Eexc plot", 50, 50, 1000, 600);
        TH2D* h_FP_Sum_HR_X_v_Eexc_coinc = histChannelPlots(reaction, "h_FP_Eex", "_HR_X_v_Eexc_coinc", c5a);
        SetTH2RangeOnCanvas(c5a, -150, 0, 0, 30);
        
        TCanvas* c5b = new TCanvas("c5b","BP Heavy residue X vs Eexc plot", 1100, 50, 1000, 600);
        TH2D* h_BP_Sum_HR_X_v_Eexc_coinc = histChannelPlots(reaction, "h_BP_Eex", "_HR_X_v_Eexc_coinc", c5b);
        SetTH2RangeOnCanvas(c5b, -150, 0, 0, 30);

        TCanvas* c5c = new TCanvas("c5c","FP Heavy residue X vs Eexc plot", 50, 700, 1000, 600);
        TH2D* hSum_HR_X_v_Eexc_coinc_theta_40_50 = histChannelPlots(reaction, "hEex", "_HR_X_v_Eexc_coinc_theta_40_50", c5c);
        SetTH2RangeOnCanvas(c5c, -150, 0, 0, 30);

        TCanvas* c5d = new TCanvas("c5d","BP Heavy residue X vs Eexc plot", 1100, 700, 1000, 600);
        TH2D* hSum_HR_X_v_Eexc_coinc_theta_10_20 = histChannelPlots(reaction, "hEex", "_HR_X_v_Eexc_coinc_theta_10_20", c5d);
        SetTH2RangeOnCanvas(c5d, -150, 0, 0, 30);
    }

    // Heavy residue coincidence fraction plots
    if(write_HR_frac) {
        TCanvas* c6a = new TCanvas("c6a","Fraction of HR transmission", 50, 50, 800, 600);
        hrFractionPlot(c6a, "h_FP_Eex");
        TCanvas* c6b = new TCanvas("c6b","Fraction of HR transmission", 900, 50, 800, 600);
        hrFractionPlot(c6b, "h_BP_Eex");
    }
    
    // Veto position plots
    if(write_veto_plots) {
        TString reaction = "206Pbdd";
        
        // Front pocket veto plots
        TCanvas* c7a = new TCanvas("c7a", "Front Pocket Veto Position Map", 50, 50, 800, 600);
        vetoPlots(c7a, "h_FP_Eex", reaction);
        
        // Back pocket veto plots
        TCanvas* c7b = new TCanvas("c7b", "Back Pocket Veto Position Map", 50, 700, 800, 600);
        vetoPlots(c7b, "h_BP_Eex", reaction);
        
        // Veto-gated banana plots
        TCanvas* c7c = new TCanvas("c7c", "Front Pocket Veto-Gated E1 Banana Plot", 900, 50, 800, 600);
        TCanvas* c7d = new TCanvas("c7d", "Front Pocket Veto-Gated Eres Banana Plot", 1750, 50, 800, 600);
        vetoBananaPlots(c7c, c7d, "h_FP_Eex", reaction);
        SetTH2RangeOnCanvas(c7c, 0, 25, 0, 10);     SetTH2RangeOnCanvas(c7d, 0, 70, 0, 30);
        
        TCanvas* c7e = new TCanvas("c7e", "Back Pocket Veto-Gated E1 Banana Plot", 900, 700, 800, 600);
        TCanvas* c7f = new TCanvas("c7f", "Back Pocket Veto-Gated Eres Banana Plot", 1750, 700, 800, 600);
        vetoBananaPlots(c7e, c7f, "h_BP_Eex", reaction);
        SetTH2RangeOnCanvas(c7e, 0, 12, 0, 2.5);     SetTH2RangeOnCanvas(c7f, 0, 100, 0, 14);
    }
}