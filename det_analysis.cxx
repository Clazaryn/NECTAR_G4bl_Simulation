// --
// NECTAR detector analysis script
// Creates TTree of processed coincidence events
// Last modification : G. Leckenby - 2025
// Based on analysis_nChamb.C and TreeAnalysis.C codes

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TMath.h>
#include <TF1.h>

#include "UtilityScripts/reaction_info.h"

// ========= Utility Functions =========

// Extract Z from PDGid
Int_t getZ(Float_t PDGid) {
    return floor((PDGid - 1000000000) * 0.0001);
}

// Extract A from PDGid
Int_t getA(Float_t PDGid) {
    Int_t Z = getZ(PDGid);
    return (PDGid - 1000000000 - Z * 10000) * 0.1;
}

// Energy resolution function
Float_t getEres(Float_t EnergyMeV, Float_t Res_Percent) {
    return gRandom->Gaus(EnergyMeV, Res_Percent / 100 * EnergyMeV);
}

// GAGG scintillator resolution
double GAGG_resolution(double E) {
    const double a = 11.2262, b = -0.573112;
    if (E <= 0) return -1.0;
    else if (E <= 1) return 10.0 / 2.355;
    else return (a * pow(E, b)) / 2.355;
}

// Read ejectile file to get true values
std::unordered_map<int, std::tuple<double, double, double, int, int>> readEjectileFile(
    const char* reaction, const char* recType, Int_t excLabel) {
    
    std::unordered_map<int, std::tuple<double, double, double, int, int>> eventMap;
    
    TString filename = Form("../%s_sim/Event_output/output_event_generator_%s_%s_excEn%02d_ejectile.txt",
                           reaction, reaction, recType, excLabel);
    
    std::ifstream file(filename.Data());
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open ejectile file: " << filename << std::endl;
        return eventMap;
    }
    
    std::string line;
    std::getline(file, line); // Skip header
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x, y, z, px, py, pz, t;
        int pdgid, eventid, trackid, parentid;
        double weight;
        
        if (iss >> x >> y >> z >> px >> py >> pz >> t >> pdgid >> eventid >> trackid >> parentid >> weight) {
            double p_mag = sqrt(px*px + py*py + pz*pz);
            double theta_true = acos(pz / p_mag) * 180.0 / M_PI;
            
            double proton_mass = 938.7830736444523; // MeV/c^2
            double E_true = sqrt(p_mag*p_mag + proton_mass*proton_mass) - proton_mass;
            double Etot_true = E_true;
            
            int Z = getZ(pdgid);
            int A = getA(pdgid);
            
            eventMap[eventid] = std::make_tuple(theta_true, E_true, Etot_true, Z, A);
        }
    }
    
    file.close();
    std::cout << "Loaded " << eventMap.size() << " events from ejectile file" << std::endl;
    return eventMap;
}

// ========= Data Classes =========

class LightEjectile {
public:
    Int_t Z, A;
    Double_t true_Estar, true_Etot, true_theta;
    Double_t recon_Estar, recon_Etot, recon_theta;
    Int_t vert_strip;
    
    LightEjectile() : Z(0), A(0), true_Estar(0), true_Etot(0), true_theta(0),
                      recon_Estar(0), recon_Etot(0), recon_theta(0), vert_strip(0) {}
};

class HeavyResidue {
public:
    Int_t Z, A;
    Double_t magsept_x, magsept_y;
    Double_t HR_x, HR_y;
    
    HeavyResidue() : Z(0), A(0), magsept_x(0), magsept_y(0), HR_x(0), HR_y(0) {}
};

class CoincEvent {
public:
    LightEjectile ejectile;
    HeavyResidue residue;
    
    CoincEvent() {}
};

// ========= Telescope Analyzer Base Class =========

class TelescopeAnalyzer {
protected:
    const char* reaction;
    const char* recType;
    Int_t excLabel;
    Double_t mass_beam, mass_target, mass_recoil, mass_ejectile;
    Double_t Ek_beam, P_beam;
    
public:
    TelescopeAnalyzer(const char* r, const char* rt, Int_t el, 
                     Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p) 
        : reaction(r), recType(rt), excLabel(el), mass_beam(mb), mass_target(mt), 
          mass_recoil(mr), mass_ejectile(me), Ek_beam(ek), P_beam(p) {}
    
    virtual ~TelescopeAnalyzer() {}
    virtual bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE, 
                              Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                              LightEjectile& ejectile) = 0;
};

// ========= New Detector Telescope Analyzer =========

class NewTelescopeAnalyzer : public TelescopeAnalyzer {
private:
    TVector3 offset;
    TVector2 rotation;
    static const Float_t Tel_DET_WIDTH = 122.0;
    static const Float_t Tel_DET_HEIGHT = 40.0;
    static const Float_t Tel_VSTRIP_WIDTH = 0.01;
    static const Float_t Tel_HSTRIP_WIDTH = 0.01;
    
public:
    NewTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                        Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                        const TVector3& off, const TVector2& rot)
        : TelescopeAnalyzer(r, rt, el, mb, mt, mr, me, ek, p), offset(off), rotation(rot) {}
    
    bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) override {
        
        // Calculate vertical strip number
        Int_t vert_strip = ceil((x_DE + Tel_DET_WIDTH/2) / Tel_VSTRIP_WIDTH);
        Int_t hor_strip = ceil((y_DE + Tel_DET_HEIGHT/2) / Tel_HSTRIP_WIDTH);
        
        if (!(1 < hor_strip && hor_strip < 40/Tel_HSTRIP_WIDTH)) {
            return false;
        }
        
        // Calculate pixel center coordinates
        Float_t x_DE_pix = -Tel_DET_WIDTH/2 + Tel_VSTRIP_WIDTH*(vert_strip-1) + Tel_VSTRIP_WIDTH/2;
        Float_t y_DE_pix = -Tel_DET_HEIGHT/2 + Tel_HSTRIP_WIDTH*(hor_strip-1) + Tel_HSTRIP_WIDTH/2;
        
        // Apply rotation to get absolute coordinates
        Float_t x_DE_abs = offset.X() + (x_DE_pix*cos(rotation.Y()) + y_DE_pix*sin(rotation.Y())*sin(rotation.X()));
        Float_t y_DE_abs = offset.Y() + (y_DE_pix*cos(rotation.X()));
        Float_t z_DE_abs = offset.Z() + (-x_DE_pix*sin(rotation.Y()) + y_DE_pix*cos(rotation.Y())*sin(rotation.X()));
        
        // Calculate theta
        TVector3 abs_vec(x_DE_abs, y_DE_abs, z_DE_abs);
        Float_t theta_calc = abs_vec.Theta() * 180.0 / M_PI;
        
        // Apply energy resolution
        Float_t Erec_DE = getEres(Edep_DE, 0.8);
        Float_t Erec_E1 = getEres(Edep_E1, 1.1);
        Float_t Erec_Eres = getEres(Edep_Eres, GAGG_resolution(Edep_Eres));
        
        // Calculate total energy
        Float_t En_tot = 0;
        if (Edep_DE != 0 && Edep_E1 == 0 && Edep_Eres == 0) {
            En_tot = Erec_DE;
        } else if (Edep_DE != 0 && Edep_E1 != 0 && Edep_Eres == 0) {
            En_tot = Erec_DE + Erec_E1;
        } else if (Edep_DE != 0 && Edep_E1 != 0 && Edep_Eres != 0) {
            En_tot = Erec_DE + Erec_E1 + Erec_Eres;
        } else {
            return false;
        }
        
        // Calculate excitation energy
        Double_t P_ejc = sqrt(pow(En_tot + mass_ejectile, 2) - pow(mass_ejectile, 2));
        Double_t Exc_en = sqrt(pow(Ek_beam + mass_beam + mass_target - En_tot - mass_ejectile, 2.) -
                               pow(P_beam, 2.) - pow(P_ejc, 2.) + 2*P_beam*P_ejc*cos(theta_calc*M_PI/180)) - mass_recoil;
        
        // Fill ejectile structure
        ejectile.vert_strip = vert_strip;
        ejectile.recon_theta = theta_calc;
        ejectile.recon_Etot = En_tot;
        ejectile.recon_Estar = Exc_en;
        
        return true;
    }
};

// ========= PoP Detector Telescope Analyzer =========

class PoPTelescopeAnalyzer : public TelescopeAnalyzer {
private:
    static const Float_t z_DE_pix = 99.025;
    TTree* tree_E2;
    std::unordered_map<int, Float_t> E2_map;
    TF1* E_corr_lowE_prot[16];
    TF1* E_corr_highE_prot[16];
    
    void initializeCorrectionFunctions() {
        // Initialize energy correction functions for each vertical strip
        // These are simplified - full implementation would include all 16 strips
        // For now, using placeholder functions
        for (int i = 0; i < 16; ++i) {
            E_corr_lowE_prot[i] = new TF1(Form("E_corr_V%02d_lowE_prot", i+1),
                "1.3+0.98*x-0.005*TMath::Power(x,2)+0.0002*TMath::Power(x,3)", 3, 18);
            E_corr_highE_prot[i] = new TF1(Form("E_corr_V%02d_highE_prot", i+1),
                "200-80*x+16*TMath::Power(x,2)-2*TMath::Power(x,3)+0.14*TMath::Power(x,4)-0.005*TMath::Power(x,5)", 3, 18);
        }
    }
    
public:
    PoPTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                         Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                         TTree* te2)
        : TelescopeAnalyzer(r, rt, el, mb, mt, mr, me, ek, p), tree_E2(te2) {
        
        initializeCorrectionFunctions();
        
        // Build E2 lookup map
        if (tree_E2) {
            Float_t EventID_E2, Edep_E2;
            tree_E2->SetBranchAddress("EventID", &EventID_E2);
            tree_E2->SetBranchAddress("Edep", &Edep_E2);
            for (Int_t i = 0; i < tree_E2->GetEntries(); ++i) {
                tree_E2->GetEntry(i);
                E2_map[static_cast<int>(EventID_E2)] = Edep_E2;
            }
        }
    }
    
    ~PoPTelescopeAnalyzer() {
        for (int i = 0; i < 16; ++i) {
            delete E_corr_lowE_prot[i];
            delete E_corr_highE_prot[i];
        }
    }
    
    bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) override {
        
        // Calculate strip numbers (PoP geometry)
        Int_t Tel_VSTRIP = 17 - ceil((x_DE + 10) / 1.25);
        Int_t Tel_HSTRIP = ceil((y_DE + 10) / 1.25);
        
        if (!(1 < Tel_HSTRIP && Tel_HSTRIP < 16) || Tel_VSTRIP < 1 || Tel_VSTRIP > 16) {
            return false;
        }
        
        // Calculate pixel coordinates
        Float_t x_DE_pix = -10 + 1.25*(Tel_VSTRIP-1) + 1.25/2;
        Float_t y_DE_pix = -10 + 1.25*(Tel_HSTRIP-1) + 1.25/2;
        
        // PoP theta calculation
        Float_t zz = z_DE_pix*cos(60*M_PI/180) + x_DE_pix*cos(30*M_PI/180);
        Float_t xx = z_DE_pix*sin(60*M_PI/180) - x_DE_pix*sin(30*M_PI/180);
        Float_t rho = sqrt(pow(xx, 2.0) + pow(y_DE_pix, 2.0) + pow(zz, 2.0));
        Float_t theta_DE_pix = acos(zz / rho);
        
        // Check if in E2
        bool isInE2 = (E2_map.count(eventID) > 0 && E2_map[eventID] > 0);
        
        // Energy reconstruction with resolution
        Float_t DSSSD = getEres(Edep_DE, 1.0);
        Float_t Etot_E = getEres(Edep_E1, 1.1);
        Float_t E_DE_E = Etot_E + DSSSD;
        
        // Apply energy correction based on strip and E2 status
        Float_t Etot_tel_res;
        if (isInE2) {
            Etot_tel_res = E_corr_highE_prot[Tel_VSTRIP-1]->Eval(E_DE_E);
        } else {
            Etot_tel_res = E_corr_lowE_prot[Tel_VSTRIP-1]->Eval(E_DE_E);
        }
        
        // Calculate excitation energy
        Double_t P_tel_res = sqrt(pow(Etot_tel_res + mass_ejectile, 2) - pow(mass_ejectile, 2));
        Double_t Exc_energy = sqrt(pow(Ek_beam + mass_beam + mass_target - Etot_tel_res - mass_ejectile, 2.) -
                                   pow(P_beam, 2.) - pow(P_tel_res, 2.) + 2*P_beam*P_tel_res*TMath::Cos(theta_DE_pix)) - mass_recoil;
        
        // Fill ejectile structure
        ejectile.vert_strip = Tel_VSTRIP;
        ejectile.recon_theta = theta_DE_pix * 180.0 / M_PI;
        ejectile.recon_Etot = Etot_tel_res;
        ejectile.recon_Estar = Exc_energy;
        
        return true;
    }
};

// ========= Main Analysis Function =========

void det_analysis(Int_t excLabel, const char* recType, const char* detType) {
    
    // Initialize reaction information
    ReactionInfo reactionInfo;
    if (!reactionInfo.loadFromIni("reac_info_nChamb.txt")) {
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
    
    if (excLabel >= reactionInfo.recoil_excEns.size()) {
        std::cerr << "Error: excLabel out of range" << std::endl;
        return;
    }
    Float_t excEn = reactionInfo.recoil_excEns[excLabel];
    
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
    TString output_filename = Form("../%s_sim/Tree_output/events_%s_%s_excEn%02d.root",
                                   reaction, reaction, recType, excLabel);
    TFile* output_file = new TFile(output_filename.Data(), "RECREATE");
    
    // Create output tree
    TTree* tree = new TTree("events", "Coincidence events");
    
    // Light ejectile branches
    Int_t ej_Z, ej_A;
    Double_t ej_true_Estar, ej_true_Etot, ej_true_theta;
    Double_t ej_recon_Estar, ej_recon_Etot, ej_recon_theta;
    Int_t ej_vert_strip;
    
    tree->Branch("ej_Z", &ej_Z);
    tree->Branch("ej_A", &ej_A);
    tree->Branch("ej_true_Estar", &ej_true_Estar);
    tree->Branch("ej_true_Etot", &ej_true_Etot);
    tree->Branch("ej_true_theta", &ej_true_theta);
    tree->Branch("ej_recon_Estar", &ej_recon_Estar);
    tree->Branch("ej_recon_Etot", &ej_recon_Etot);
    tree->Branch("ej_recon_theta", &ej_recon_theta);
    tree->Branch("ej_vert_strip", &ej_vert_strip);
    
    // Heavy residue branches
    Int_t hr_Z, hr_A;
    Double_t hr_magsept_x, hr_magsept_y;
    Double_t hr_HR_x, hr_HR_y;
    
    tree->Branch("hr_Z", &hr_Z);
    tree->Branch("hr_A", &hr_A);
    tree->Branch("hr_magsept_x", &hr_magsept_x);
    tree->Branch("hr_magsept_y", &hr_magsept_y);
    tree->Branch("hr_HR_x", &hr_HR_x);
    tree->Branch("hr_HR_y", &hr_HR_y);
    
    // Get detector trees based on detector type
    TelescopeAnalyzer* analyzer = nullptr;
    
    if (strcmp(detType, "New") == 0) {
        // New detector setup
        TTree* tree_DE = (TTree*)ejectile_file->Get("Detector/Front_tele_DE");
        TTree* tree_E1 = (TTree*)ejectile_file->Get("Detector/Front_tele_E1");
        TTree* tree_Eres = (TTree*)ejectile_file->Get("Detector/Front_tele_Eres");
        
        if (!tree_DE || !tree_E1 || !tree_Eres) {
            std::cerr << "Error: Could not find New detector trees" << std::endl;
            return;
        }
        
        TVector3 offset(-80.858, 49.802, 31.162);
        TVector2 rotation(-20*M_PI/180, -20*M_PI/180);
        analyzer = new NewTelescopeAnalyzer(reaction, recType, excLabel,
                                           mass_beam, mass_target, mass_recoil, mass_ejectile,
                                           Ek_beam, P_beam, offset, rotation);
        
        // Build lookup maps
        std::unordered_map<int, Float_t> E1_map, Eres_map;
        Float_t EventID_E1, Edep_E1, EventID_Eres, Edep_Eres;
        tree_E1->SetBranchAddress("EventID", &EventID_E1);
        tree_E1->SetBranchAddress("Edep", &Edep_E1);
        tree_Eres->SetBranchAddress("EventID", &EventID_Eres);
        tree_Eres->SetBranchAddress("Edep", &Edep_Eres);
        
        for (Int_t i = 0; i < tree_E1->GetEntries(); ++i) {
            tree_E1->GetEntry(i);
            E1_map[static_cast<int>(EventID_E1)] = Edep_E1;
        }
        for (Int_t i = 0; i < tree_Eres->GetEntries(); ++i) {
            tree_Eres->GetEntry(i);
            Eres_map[static_cast<int>(EventID_Eres)] = Edep_Eres;
        }
        
        // Get HR trees
        TTree* dipole_testblock = (TTree*)recoil_file->Get("VirtualDetector/dipole_testblock");
        TTree* HR_testblock = (TTree*)recoil_file->Get("Detector/HR_testblock");
        
        std::unordered_map<int, std::tuple<Float_t, Float_t, Float_t>> Dip_map, HR_map;
        Float_t EventID_Dip, x_Dip, y_Dip, z_Dip;
        Float_t EventID_HR, x_HR, y_HR, z_HR;
        
        if (dipole_testblock) {
            dipole_testblock->SetBranchAddress("EventID", &EventID_Dip);
            dipole_testblock->SetBranchAddress("x", &x_Dip);
            dipole_testblock->SetBranchAddress("y", &y_Dip);
            dipole_testblock->SetBranchAddress("z", &z_Dip);
            for (Int_t i = 0; i < dipole_testblock->GetEntries(); ++i) {
                dipole_testblock->GetEntry(i);
                Dip_map[static_cast<int>(EventID_Dip)] = std::make_tuple(x_Dip, y_Dip, z_Dip);
            }
        }
        
        std::unordered_map<int, Float_t> HR_PDGid_map;
        Float_t PDGid_HR;
        if (HR_testblock) {
            HR_testblock->SetBranchAddress("EventID", &EventID_HR);
            HR_testblock->SetBranchAddress("x", &x_HR);
            HR_testblock->SetBranchAddress("y", &y_HR);
            HR_testblock->SetBranchAddress("z", &z_HR);
            if (HR_testblock->GetBranch("PDGid")) {
                HR_testblock->SetBranchAddress("PDGid", &PDGid_HR);
            }
            for (Int_t i = 0; i < HR_testblock->GetEntries(); ++i) {
                HR_testblock->GetEntry(i);
                HR_map[static_cast<int>(EventID_HR)] = std::make_tuple(x_HR, y_HR, z_HR);
                if (HR_testblock->GetBranch("PDGid")) {
                    HR_PDGid_map[static_cast<int>(EventID_HR)] = PDGid_HR;
                }
            }
        }
        
        // Process events
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
            
            LightEjectile ejectile;
            if (analyzer->analyzeEvent(eventID, x_DE, y_DE, z_DE, Edep_DE, Edep_E1_val, Edep_Eres_val, ejectile)) {
                // Fill true values if available
                if (eventMap.count(eventID) > 0) {
                    auto true_vals = eventMap[eventID];
                    ejectile.true_theta = std::get<0>(true_vals);
                    ejectile.true_Etot = std::get<2>(true_vals);
                    ejectile.true_Estar = std::get<1>(true_vals); // Using E_true as Estar approximation
                    ejectile.Z = std::get<3>(true_vals);
                    ejectile.A = std::get<4>(true_vals);
                }
                
                // Check for HR coincidence
                if (HR_map.count(eventID) > 0 || Dip_map.count(eventID) > 0) {
                    HeavyResidue residue;
                    
                    if (Dip_map.count(eventID) > 0) {
                        auto dip_pos = Dip_map[eventID];
                        residue.magsept_x = std::get<0>(dip_pos);
                        residue.magsept_y = std::get<1>(dip_pos);
                    }
                    
                    if (HR_map.count(eventID) > 0) {
                        auto hr_pos = HR_map[eventID];
                        residue.HR_x = std::get<0>(hr_pos);
                        residue.HR_y = std::get<1>(hr_pos);
                        
                        // Extract Z and A from PDGid if available
                        if (HR_PDGid_map.count(eventID) > 0) {
                            Float_t pdgid = HR_PDGid_map[eventID];
                            residue.Z = getZ(pdgid);
                            residue.A = getA(pdgid);
                        }
                    }
                    
                    // Fill tree
                    ej_Z = ejectile.Z;
                    ej_A = ejectile.A;
                    ej_true_Estar = ejectile.true_Estar;
                    ej_true_Etot = ejectile.true_Etot;
                    ej_true_theta = ejectile.true_theta;
                    ej_recon_Estar = ejectile.recon_Estar;
                    ej_recon_Etot = ejectile.recon_Etot;
                    ej_recon_theta = ejectile.recon_theta;
                    ej_vert_strip = ejectile.vert_strip;
                    
                    hr_Z = residue.Z;
                    hr_A = residue.A;
                    hr_magsept_x = residue.magsept_x;
                    hr_magsept_y = residue.magsept_y;
                    hr_HR_x = residue.HR_x;
                    hr_HR_y = residue.HR_y;
                    
                    tree->Fill();
                }
            }
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
        
        analyzer = new PoPTelescopeAnalyzer(reaction, recType, excLabel,
                                           mass_beam, mass_target, mass_recoil, mass_ejectile,
                                           Ek_beam, P_beam, tree_E2);
        
        // Build lookup maps
        std::unordered_map<int, Float_t> E1_map;
        Float_t EventID_E1, Edep_E1;
        tree_E1->SetBranchAddress("EventID", &EventID_E1);
        tree_E1->SetBranchAddress("Edep", &Edep_E1);
        
        for (Int_t i = 0; i < tree_E1->GetEntries(); ++i) {
            tree_E1->GetEntry(i);
            E1_map[static_cast<int>(EventID_E1)] = Edep_E1;
        }
        
        // Get HR trees (PoP naming may differ)
        TTree* Virtual_Magsept = (TTree*)recoil_file->Get("VirtualDetector/MagSept_virtual");
        TTree* Virtual_HR = (TTree*)recoil_file->Get("VirtualDetector/Virtual_3n");
        
        std::unordered_map<int, std::tuple<Float_t, Float_t, Float_t>> Magsept_map, HR_map;
        Float_t EventID_magsept, x_magsept, y_magsept, z_magsept;
        Float_t EventID_virtual, x_virtual, y_virtual, z_virtual;
        
        if (Virtual_Magsept) {
            Virtual_Magsept->SetBranchAddress("EventID", &EventID_magsept);
            Virtual_Magsept->SetBranchAddress("x", &x_magsept);
            Virtual_Magsept->SetBranchAddress("y", &y_magsept);
            Virtual_Magsept->SetBranchAddress("z", &z_magsept);
            for (Int_t i = 0; i < Virtual_Magsept->GetEntries(); ++i) {
                Virtual_Magsept->GetEntry(i);
                Magsept_map[static_cast<int>(EventID_magsept)] = std::make_tuple(x_magsept, y_magsept, z_magsept);
            }
        }
        
        std::unordered_map<int, Float_t> HR_PDGid_map;
        Float_t PDGid_virtual;
        if (Virtual_HR) {
            Virtual_HR->SetBranchAddress("EventID", &EventID_virtual);
            Virtual_HR->SetBranchAddress("x", &x_virtual);
            Virtual_HR->SetBranchAddress("y", &y_virtual);
            Virtual_HR->SetBranchAddress("z", &z_virtual);
            if (Virtual_HR->GetBranch("PDGid")) {
                Virtual_HR->SetBranchAddress("PDGid", &PDGid_virtual);
            }
            for (Int_t i = 0; i < Virtual_HR->GetEntries(); ++i) {
                Virtual_HR->GetEntry(i);
                HR_map[static_cast<int>(EventID_virtual)] = std::make_tuple(x_virtual, y_virtual, z_virtual);
                if (Virtual_HR->GetBranch("PDGid")) {
                    HR_PDGid_map[static_cast<int>(EventID_virtual)] = PDGid_virtual;
                }
            }
        }
        
        // Process events
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
            
            LightEjectile ejectile;
            if (analyzer->analyzeEvent(eventID, x_DE, y_DE, z_DE, Edep_DE, Edep_E1_val, 0, ejectile)) {
                // Fill true values if available
                if (eventMap.count(eventID) > 0) {
                    auto true_vals = eventMap[eventID];
                    ejectile.true_theta = std::get<0>(true_vals);
                    ejectile.true_Etot = std::get<2>(true_vals);
                    ejectile.true_Estar = std::get<1>(true_vals);
                    ejectile.Z = std::get<3>(true_vals);
                    ejectile.A = std::get<4>(true_vals);
                }
                
                // Check for HR coincidence
                if (HR_map.count(eventID) > 0 || Magsept_map.count(eventID) > 0) {
                    HeavyResidue residue;
                    
                    if (Magsept_map.count(eventID) > 0) {
                        auto mag_pos = Magsept_map[eventID];
                        residue.magsept_x = std::get<0>(mag_pos);
                        residue.magsept_y = std::get<1>(mag_pos);
                    }
                    
                    if (HR_map.count(eventID) > 0) {
                        auto hr_pos = HR_map[eventID];
                        residue.HR_x = std::get<0>(hr_pos);
                        residue.HR_y = std::get<1>(hr_pos);
                        
                        // Extract Z and A from PDGid if available
                        if (HR_PDGid_map.count(eventID) > 0) {
                            Float_t pdgid = HR_PDGid_map[eventID];
                            residue.Z = getZ(pdgid);
                            residue.A = getA(pdgid);
                        }
                    }
                    
                    // Fill tree
                    ej_Z = ejectile.Z;
                    ej_A = ejectile.A;
                    ej_true_Estar = ejectile.true_Estar;
                    ej_true_Etot = ejectile.true_Etot;
                    ej_true_theta = ejectile.true_theta;
                    ej_recon_Estar = ejectile.recon_Estar;
                    ej_recon_Etot = ejectile.recon_Etot;
                    ej_recon_theta = ejectile.recon_theta;
                    ej_vert_strip = ejectile.vert_strip;
                    
                    hr_Z = residue.Z;
                    hr_A = residue.A;
                    hr_magsept_x = residue.magsept_x;
                    hr_magsept_y = residue.magsept_y;
                    hr_HR_x = residue.HR_x;
                    hr_HR_y = residue.HR_y;
                    
                    tree->Fill();
                }
            }
        }
    } else {
        std::cerr << "Error: Unknown detector type: " << detType << std::endl;
        return;
    }
    
    // Write and close
    output_file->cd();
    tree->Write();
    output_file->Close();
    
    delete analyzer;
    ejectile_file->Close();
    recoil_file->Close();
    
    std::cout << "Analysis complete. Output written to: " << output_filename << std::endl;
}
