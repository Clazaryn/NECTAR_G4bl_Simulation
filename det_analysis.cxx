// --
// NECTAR detector analysis script
// Creates TTree of processed coincidence events
// Last modification : G. Leckenby - 2025
// Based on analysis_nChamb.C and TreeAnalysis.C codes

#include "UtilityScripts/det_analysis.h"
#include "UtilityScripts/reaction_info.h"

// ========= Helper Function to Process Telescope Events =========

void processTelescopeEvents(TTree* tree_DE, TTree* tree_E1, TTree* tree_Eres,
                            TelescopeAnalyzer* analyzer,
                            const std::unordered_map<int, std::tuple<double, double, double, int, int>>& eventMap,
                            const VirtualDetectorData& magsept_data,
                            const VirtualDetectorData& hrplane_data,
                            const VirtualDetectorData& quadwall_data,
                            Int_t detector_id,
                            TTree* output_tree,
                            LightEjectile* ejectile,
                            HeavyResidue* residue) {
    
    // Build lookup maps for E detectors
    std::unordered_map<int, Float_t> E1_map, Eres_map;
    Float_t EventID_E1, Edep_E1, EventID_Eres, Edep_Eres;
    
    if (tree_E1) {
        tree_E1->SetBranchAddress("EventID", &EventID_E1);
        tree_E1->SetBranchAddress("Edep", &Edep_E1);
        for (Int_t i = 0; i < tree_E1->GetEntries(); ++i) {
            tree_E1->GetEntry(i);
            E1_map[static_cast<int>(EventID_E1)] = Edep_E1;
        }
    }
    
    if (tree_Eres) {
        tree_Eres->SetBranchAddress("EventID", &EventID_Eres);
        tree_Eres->SetBranchAddress("Edep", &Edep_Eres);
        for (Int_t i = 0; i < tree_Eres->GetEntries(); ++i) {
            tree_Eres->GetEntry(i);
            Eres_map[static_cast<int>(EventID_Eres)] = Edep_Eres;
        }
    }
    
    // Process DE events
    Float_t x_DE, y_DE, z_DE, EventID_DE, Edep_DE;
    tree_DE->SetBranchAddress("x", &x_DE);
    tree_DE->SetBranchAddress("y", &y_DE);
    tree_DE->SetBranchAddress("z", &z_DE);
    tree_DE->SetBranchAddress("EventID", &EventID_DE);
    tree_DE->SetBranchAddress("Edep", &Edep_DE);
    
    for (Int_t i = 0; i < tree_DE->GetEntries(); ++i) {
        tree_DE->GetEntry(i);
        Int_t eventID = static_cast<int>(EventID_DE);
        
        Float_t Edep_E1_val = (E1_map.count(eventID) > 0) ? E1_map[eventID] : 0;
        Float_t Edep_Eres_val = (Eres_map.count(eventID) > 0) ? Eres_map[eventID] : 0;
        
        if (analyzer->analyzeEvent(eventID, x_DE, y_DE, z_DE, Edep_DE, Edep_E1_val, Edep_Eres_val, *ejectile)) {
            // Fill true values if available
            auto it = eventMap.find(eventID);
            if (it != eventMap.end()) {
                auto true_vals = it->second;
                ejectile->true_theta = std::get<0>(true_vals);
                ejectile->true_Etot = std::get<2>(true_vals);
                ejectile->true_Estar = std::get<1>(true_vals);
                ejectile->Z = std::get<3>(true_vals);
                ejectile->A = std::get<4>(true_vals);
            }
            
            // Set detector ID
            ejectile->detector_id = detector_id;
            
            // Check for HR coincidence (need at least MagSept)
            if (magsept_data.pos_map.count(eventID) > 0) {
                fillResidueFromMaps(eventID, magsept_data, hrplane_data, quadwall_data, residue);
                output_tree->Fill();
            }
        }
    }
}

// ========= Main Analysis Function =========

void det_analysis(Int_t excLabel, const char* recType) {
    
    // Initialize reaction information
    ReactionInfo reactionInfo;
    if (!reactionInfo.loadFromIni("reac_info.txt")) {
        std::cerr << "Error: Could not load reaction info. Exiting." << std::endl;
        return;
    }
    
    const char* reaction = reactionInfo.reaction.c_str();
    Double_t mass_beam = reactionInfo.mass_beam;
    Double_t mass_target = reactionInfo.mass_target;
    Double_t mass_recoil = reactionInfo.mass_recoil;
    Double_t mass_ejectile = reactionInfo.mass_eject;
    Float_t EuA_beam = reactionInfo.beam_EuA;
    Double_t Ek_beam = mass_beam / 931.49410242 * EuA_beam;
    Double_t P_beam = sqrt(pow(Ek_beam + mass_beam, 2.0) - pow(mass_beam, 2.0));
    
    // Determine detector type from INI file
    std::string det_setup = reactionInfo.det_setup;
    const char* detType;
    if (det_setup == "new" || det_setup == "New") {
        detType = "New";
    } else if (det_setup == "PoP" || det_setup == "pop") {
        detType = "PoP";
    } else {
        std::cerr << "Error: Unknown detector setup in INI file: " << det_setup << std::endl;
        std::cerr << "Expected 'new' or 'PoP'" << std::endl;
        return;
    }
    
    if (excLabel >= reactionInfo.recoil_excEns.size()) {
        std::cerr << "Error: excLabel out of range" << std::endl;
        return;
    }
    
    // Read true values from ejectile file
    auto eventMap = readEjectileFile(reaction, recType, excLabel);
    
    // Open input files
    TString ejectile_filename = Form("../%s_sim/Detector_output/Detectors_%s_%s_excEn%02d_ejectile.root",
                                    reaction, reaction, recType, excLabel);
    TString recoil_filename = Form("../%s_sim/Detector_output/Detectors_%s_%s_excEn%02d_recoil.root",
                                   reaction, reaction, recType, excLabel);
    
    TFile* ejectile_file = TFile::Open(ejectile_filename.Data());
    TFile* recoil_file = TFile::Open(recoil_filename.Data());
    
    if (!ejectile_file || !recoil_file) {
        std::cerr << "Error: Could not open input files" << std::endl;
        return;
    }
    
    // Create output file
    TString output_filename = Form("../%s_sim/Det_analysis/events_%s_%s_excEn%02d.root",
                                   reaction, reaction, recType, excLabel);
    TFile* output_file = new TFile(output_filename.Data(), "RECREATE");
    
    // Create output tree
    TTree* tree = new TTree("events", "Coincidence events");
    
    // Create objects for branches
    LightEjectile* ejectile = new LightEjectile();
    HeavyResidue* residue = new HeavyResidue();
    
    // Set branches using objects
    tree->Branch("ejectile", &ejectile);
    tree->Branch("residue", &residue);
    
    // Load virtual detector data (same names for both detector types)
    VirtualDetectorData magsept_data = loadVirtualDetector(recoil_file, "VirtualDetector/virt_MagSept");
    VirtualDetectorData hrplane_data = loadVirtualDetector(recoil_file, "VirtualDetector/virt_HRplane");
    VirtualDetectorData quadwall_data = loadVirtualDetector(recoil_file, "VirtualDetector/virt_QuadWall");
    
    if (strcmp(detType, "New") != 0 && strcmp(detType, "PoP") != 0) {
        std::cerr << "Error: Unknown detector type: " << detType << std::endl;
        return;
    }
    
    // Process detector-specific telescopes
    if (strcmp(detType, "New") == 0) {
        // New detector setup - process both primary and auxillary telescopes
        std::vector<TString> detector_names = {"primary", "auxillary"};
        std::vector<TVector3> offsets = {
            TVector3(-80.858, 49.802, 31.162),      // primary
            TVector3(0, 49.802, 206.54)             // auxillary
        };
        std::vector<TVector2> rotations = {
            TVector2(-20*M_PI/180, -20*M_PI/180),   // primary
            TVector2(-20*M_PI/180, 0)               // auxillary
        };
        
        // Process each telescope
        for (size_t tel_idx = 0; tel_idx < detector_names.size(); ++tel_idx) {
            TString tel_name = detector_names[tel_idx];
            
            TTree* tree_DE = (TTree*)ejectile_file->Get(Form("Detector/%s_tele_DE", tel_name.Data()));
            TTree* tree_E1 = (TTree*)ejectile_file->Get(Form("Detector/%s_tele_E1", tel_name.Data()));
            TTree* tree_Eres = (TTree*)ejectile_file->Get(Form("Detector/%s_tele_Eres", tel_name.Data()));
            
            if (!tree_DE || !tree_E1 || !tree_Eres) {
                std::cerr << "Warning: Could not find " << tel_name << " detector trees" << std::endl;
                continue;
            }
            
            NewTelescopeAnalyzer analyzer(reaction, recType, excLabel,
                                          mass_beam, mass_target, mass_recoil, mass_ejectile,
                                          Ek_beam, P_beam, offsets[tel_idx], rotations[tel_idx]);
            
            processTelescopeEvents(tree_DE, tree_E1, tree_Eres, &analyzer, eventMap,
                                  magsept_data, hrplane_data, quadwall_data, tel_idx, tree, ejectile, residue);
        }
        
    } else if (strcmp(detType, "PoP") == 0) {
        // PoP detector setup
        TTree* tree_DE = (TTree*)ejectile_file->Get("Detector/Telescope_DE");
        TTree* tree_E1 = (TTree*)ejectile_file->Get("Detector/Telescope_E1");
        TTree* tree_E2 = (TTree*)ejectile_file->Get("Detector/Telescope_E2");
        
        if (!tree_DE || !tree_E1 || !tree_E2) {
            std::cerr << "Error: Could not find PoP detector trees" << std::endl;
            return;
        }
        
        PoPTelescopeAnalyzer analyzer(reaction, recType, excLabel,
                                     mass_beam, mass_target, mass_recoil, mass_ejectile,
                                     Ek_beam, P_beam, tree_E2);
        
        // For PoP, Eres is not used (set to nullptr)
        processTelescopeEvents(tree_DE, tree_E1, nullptr, &analyzer, eventMap,
                              magsept_data, hrplane_data, quadwall_data, 0, tree, ejectile, residue);
    }
    
    // Write and close
    output_file->cd();
    tree->Write();
    output_file->Close();
    
    ejectile_file->Close();
    recoil_file->Close();
    
    delete ejectile;
    delete residue;
    
    std::cout << "Analysis complete. Output written to: " << output_filename << std::endl;
}
