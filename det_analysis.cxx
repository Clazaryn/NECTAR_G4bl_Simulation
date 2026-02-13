// --
// NECTAR simulation: detector analysis script
// Creates TTree of processed coincidence events
// Last modification : G. Leckenby - 2025
// Based on analysis_nChamb.C and TreeAnalysis.C codes

#include "UtilityScripts/det_analysis.h"
#include <cstdlib>

// ========= Helper Function to Process Telescope Events =========
// this function implements the analysis of the telescope events into reconstructed physics values

void processTelescopeEvents(TTree* tree_DE, TTree* tree_E1, TTree* tree_Eres,
                            TelescopeAnalyzer* analyzer,
                            const std::unordered_map<int, std::tuple<double, double, double, int, int>>& eventMap,
                            const VirtualDetectorData& magsept_data,
                            const VirtualDetectorData& hrplane_data,
                            const VirtualDetectorData& quadwall_data,
                            Int_t detector_id,
                            TTree* output_tree,
                            LightEjectile* ejectile,
                            HeavyResidue* residue,
                            TFile* eject_det_output = nullptr) {
    
    // Build lookup maps for E detectors - these avoid recursive for loops
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
    
    // For PoP: individual E detector maps (E2 through E6)
    std::unordered_map<int, Float_t> E2_map, E3_map, E4_map, E5_map, E6_map;
    
    if (tree_Eres) {
        // New detector setup: single Eres detector (CsI(Tl) crystal) called Eres, so we create lookup map directly
        tree_Eres->SetBranchAddress("EventID", &EventID_Eres);
        tree_Eres->SetBranchAddress("Edep", &Edep_Eres);
        for (Int_t i = 0; i < tree_Eres->GetEntries(); ++i) {
            tree_Eres->GetEntry(i);
            Eres_map[static_cast<int>(EventID_Eres)] = Edep_Eres;
        }
    } else if (eject_det_output) {
        // PoP detector setup: store individual E2-E6 values in separate maps
        // Also sum them for backward compatibility with analyzer (which uses sum for reconstruction)
        Float_t EventID_E, Edep_E;
        std::unordered_map<int, Float_t>* E_maps[] = {&E2_map, &E3_map, &E4_map, &E5_map, &E6_map};
        for (int i = 2; i <= 6; ++i) {  // E2 through E6
            TString tree_name = Form("Detector/Telescope_E%d", i);
            TTree* tree_E_i = (TTree*)eject_det_output->Get(tree_name.Data());
            if (tree_E_i) {
                tree_E_i->SetBranchAddress("EventID", &EventID_E);
                tree_E_i->SetBranchAddress("Edep", &Edep_E);
                std::unordered_map<int, Float_t>* current_map = E_maps[i-2];  // i-2 because array starts at E2
                for (Int_t j = 0; j < tree_E_i->GetEntries(); ++j) {
                    tree_E_i->GetEntry(j);
                    int eventID = static_cast<int>(EventID_E);
                    (*current_map)[eventID] = Edep_E;  // Store individual value
                    Eres_map[eventID] += Edep_E;      // Also sum for analyzer compatibility
                }
            }
        }
    }
    
    // Process DE events from DSSD detector
    Float_t x_DE, y_DE, z_DE, EventID_DE, Edep_DE;
    tree_DE->SetBranchAddress("x", &x_DE);
    tree_DE->SetBranchAddress("y", &y_DE);
    tree_DE->SetBranchAddress("z", &z_DE);
    tree_DE->SetBranchAddress("EventID", &EventID_DE);
    tree_DE->SetBranchAddress("Edep", &Edep_DE);
    
    for (Int_t i = 0; i < tree_DE->GetEntries(); ++i) {
        tree_DE->GetEntry(i);
        Int_t eventID = static_cast<int>(EventID_DE);
        
        // Build vector of deposited energies: [dE, E1, E2, E3, E4, E5, E6]
        // For New detectors: E2-E6 will be 0, and Eres is passed separately (handled in analyzeEvent)
        // For PoP detectors: all values are populated
        std::vector<Float_t> Edep_vec;
        Edep_vec.push_back(Edep_DE);  // Index 0: dE
        Edep_vec.push_back((E1_map.count(eventID) > 0) ? E1_map[eventID] : 0);  // Index 1: E1
        
        // For PoP detectors: loop through E_maps to populate E2-E6
        if (!tree_Eres && eject_det_output) {
            std::unordered_map<int, Float_t>* E_maps[] = {&E2_map, &E3_map, &E4_map, &E5_map, &E6_map};
            for (int j = 0; j < 5; ++j) {  // E2 through E6
                Edep_vec.push_back((E_maps[j]->count(eventID) > 0) ? (*E_maps[j])[eventID] : 0);
            }
        } else {
            // For New detectors: E2-E6 are 0, but we still need to pass Eres
            // We'll pass it as E2 (index 2), and E3-E6 will be 0
            Float_t Edep_Eres_val = (Eres_map.count(eventID) > 0) ? Eres_map[eventID] : 0;
            Edep_vec.push_back(Edep_Eres_val);  // Index 2: Eres (for New detectors)
            for (int j = 0; j < 4; ++j) {  // E3-E6 are 0 for New detectors
                Edep_vec.push_back(0);
            }
        }
        
        // analyseEvent function performs the reconstruction analysis
        // tree is filled if event is successfully reconstructed and if event is hit the MagSept detector
        if (analyzer->analyzeEvent(eventID, x_DE, y_DE, z_DE, Edep_vec, *ejectile)) {
            // Fill true values if map provided from ejectile event file (Event_output)
            auto it = eventMap.find(eventID);
            if (it != eventMap.end()) {
                auto true_vals = it->second;
                ejectile->Z = std::get<0>(true_vals);
                ejectile->A = std::get<1>(true_vals);
                ejectile->true_Eexc = std::get<2>(true_vals);
                ejectile->true_Eejc = std::get<3>(true_vals);
                ejectile->true_theta = std::get<4>(true_vals);
            }
            
            ejectile->detector_id = detector_id;        // Set detector ID
            
            // Check for HR coincidence (need at least MagSept)
            if (magsept_data.pos_map.count(eventID) > 0) {
                fillResidueFromMaps(eventID, magsept_data, hrplane_data, quadwall_data, residue);
                output_tree->Fill();        // Fill tree if event is hit the MagSept detector
            }
        }
    }
}

// ========= Main Analysis Function =========

void det_analysis(Int_t excLabel, const char* recType) {
    
    // Initialize reaction information
    ReactionInfo reactionInfo;
    if (!reactionInfo.loadFromIni("reac_info.txt")) {
        std::cerr << "det_analysis.cxx:" << __LINE__ << ": Error: Could not load reaction info. Exiting." << std::endl;
        std::exit(1);
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
        std::cerr << "det_analysis.cxx:" << __LINE__ << ": Error: Unknown detector setup in INI file: " << det_setup << std::endl;
        std::cerr << "Expected 'new' or 'PoP'" << std::endl;
        std::exit(1);
    }
    
    if (excLabel >= reactionInfo.recoil_excEns.size()) {
        std::cerr << "det_analysis.cxx:" << __LINE__ << ": Error: excLabel out of range" << std::endl;
        std::exit(1);
    }
    
    // Get excitation energy value to format label (matches bash script: excEnXX.XMeV)
    Double_t excEn = reactionInfo.recoil_excEns[excLabel];
    TString lbl = Form("%04.1fMeV", excEn);  // Format: XX.XMeV (1 decimal place with leading zero)
    
    // Read true values from ejectile event file (Event_output text file)
    auto eventMap = readEjectileFile(reaction, recType, excLabel, excEn, reactionInfo);
    
    // Open input files from G4bl (format matches bash script: excEnXX.XMeV)
    TString eject_det_output_filename = Form("../%s_sim/Detector_output/Detectors_%s_%s_excEn%s_ejectile.root",
                                             reaction, reaction, recType, lbl.Data());
    TString recoil_det_output_filename = Form("../%s_sim/Detector_output/Detectors_%s_%s_excEn%s_recoil.root",
                                                reaction, reaction, recType, lbl.Data());
    
    TFile* eject_det_output = TFile::Open(eject_det_output_filename.Data());
    TFile* recoil_det_output = TFile::Open(recoil_det_output_filename.Data());
    
    if (!eject_det_output || !recoil_det_output) {
        std::cerr << "det_analysis.cxx:" << __LINE__ << ": Error: Could not open input files" << std::endl;
        std::exit(1);
    }
    
    // Create output file (use same label format for consistency)
    TString output_filename = Form("../%s_sim/Det_analysis/events_%s_%s_excEn%s.root",
                                   reaction, reaction, recType, lbl.Data());
    TFile* output_file = new TFile(output_filename.Data(), "RECREATE");
    
    // Create output tree
    TTree* tree = new TTree("events", "Coincidence events");
    
    // Create Short_t for decay channel encoding
    // Encoding: 0=HRg, 1=HR1n, 2=HR2n, 3=HR3n, 4=HR4n, -1=unknown
    Short_t decay_channel_id = -1;
    if (strcmp(recType, "HRg") == 0) {
        decay_channel_id = 0;
    } else if (strcmp(recType, "HR1n") == 0) {
        decay_channel_id = 1;
    } else if (strcmp(recType, "HR2n") == 0) {
        decay_channel_id = 2;
    } else if (strcmp(recType, "HR3n") == 0) {
        decay_channel_id = 3;
    } else if (strcmp(recType, "HR4n") == 0) {
        decay_channel_id = 4;
    }
    
    // Create objects for branches
    LightEjectile* ejectile = new LightEjectile();
    HeavyResidue* residue = new HeavyResidue();
    
    // Set branches using objects
    tree->Branch("ejectile", &ejectile);
    tree->Branch("residue", &residue);
    tree->Branch("decay_channel", &decay_channel_id);
    
    // Load virtual detector data (same names for both detector types)
    // Read PDGid from Event_output text file to avoid precision loss from ROOT Float_t storage
    VirtualDetectorData magsept_data = loadVirtualDetector(recoil_det_output, "VirtualDetector/virt_MagSept",
                                                            reaction, recType, excLabel, excEn);
    VirtualDetectorData hrplane_data = loadVirtualDetector(recoil_det_output, "VirtualDetector/virt_HRplane",
                                                            reaction, recType, excLabel, excEn);
    VirtualDetectorData quadwall_data = loadVirtualDetector(recoil_det_output, "VirtualDetector/virt_QuadWall",
                                                            reaction, recType, excLabel, excEn);
    
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
            
            TTree* tree_DE = (TTree*)eject_det_output->Get(Form("Detector/%s_tele_DE", tel_name.Data()));
            TTree* tree_E1 = (TTree*)eject_det_output->Get(Form("Detector/%s_tele_E1", tel_name.Data()));
            TTree* tree_Eres = (TTree*)eject_det_output->Get(Form("Detector/%s_tele_Eres", tel_name.Data()));
            
            if (!tree_DE || !tree_E1 || !tree_Eres) {
                std::cerr << "Warning: Could not find " << tel_name << " detector trees" << std::endl;
                continue;
            }
            
            NewTelescopeAnalyzer analyzer(reaction, recType, excLabel,
                                          mass_beam, mass_target, mass_recoil, mass_ejectile,
                                          Ek_beam, P_beam, offsets[tel_idx], rotations[tel_idx]);
            
            processTelescopeEvents(tree_DE, tree_E1, tree_Eres, &analyzer, eventMap,
                                  magsept_data, hrplane_data, quadwall_data, tel_idx, tree, ejectile, residue, eject_det_output);
        }
        
    } else if (strcmp(detType, "PoP") == 0) {
        // PoP detector setup: positions, etc. are hardcoded
        TTree* tree_DE = (TTree*)eject_det_output->Get("Detector/Telescope_DE");
        TTree* tree_E1 = (TTree*)eject_det_output->Get("Detector/Telescope_E1");
        TTree* tree_E2 = (TTree*)eject_det_output->Get("Detector/Telescope_E2");
        
        if (!tree_DE || !tree_E1 || !tree_E2) {
            std::cerr << "det_analysis.cxx:" << __LINE__ << ": Error: Could not find PoP detector trees" << std::endl;
            std::exit(1);
        }
        
        PoPTelescopeAnalyzer analyzer(reaction, recType, excLabel,
                                     mass_beam, mass_target, mass_recoil, mass_ejectile,
                                     Ek_beam, P_beam, tree_E2);
        
        // For PoP, Eres is sum of E2+E3+E4+... (calculated in processTelescopeEvents)
        processTelescopeEvents(tree_DE, tree_E1, nullptr, &analyzer, eventMap,
                              magsept_data, hrplane_data, quadwall_data, 0, tree, ejectile, residue, eject_det_output);
    }
    
    // Write and close
    output_file->cd();
    tree->Write();
    output_file->Close();
    
    eject_det_output->Close();
    recoil_det_output->Close();
    
    delete ejectile;
    delete residue;
    
    std::cout << "Analysis complete. Output written to: " << output_filename << std::endl;
}
