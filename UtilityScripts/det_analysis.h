#ifndef DET_ANALYSIS_H
#define DET_ANALYSIS_H

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
#include <TVector3.h>
#include <TVector2.h>
#include <TMath.h>
#include <TF1.h>

#include "reaction_info.h"

// ========= Data Classes =========

class LightEjectile {
public:
    Int_t Z, A;
    Int_t detector_id;  // 0=primary, 1=auxillary for New detectors
    Int_t vert_strip, hor_strip;
    Double_t true_Eexc, true_Eejc, true_theta;
    Double_t recon_Eexc, recon_Eejc, recon_theta;
    Float_t meas_dE, meas_E1, meas_Eres;  // Measured Eres energy (crystal in New setup)
    Float_t meas_E2, meas_E3, meas_E4, meas_E5, meas_E6;  // Measured E2 through E6 energy (thick Si in PoP setup)
    
    LightEjectile();
};

class HeavyResidue {
public:
    Int_t Z, A;
    Bool_t hit_MagSept;  // true if event hit MagSept virtual detector
    Double_t MagSept_x, MagSept_y;  // x,y coordinates of event in MagSept virtual detector
    Bool_t hit_HRplane;  // true if event hit HRplane virtual detector
    Double_t HRplane_x, HRplane_y;  // x,y coordinates of event in HRplane virtual detector
    Bool_t hit_QuadWall;  // true if event hit QuadWall virtual detector
    Double_t QuadWall_x, QuadWall_y;  // x,y coordinates of event in QuadWall virtual detector
    
    HeavyResidue();
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
                     Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p);
    virtual ~TelescopeAnalyzer() {}
    virtual bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE, 
                              const std::vector<Float_t>& Edep_vec,
                              LightEjectile& ejectile) = 0;
};

// ========= New Detector Telescope Analyzer =========

class NewTelescopeAnalyzer : public TelescopeAnalyzer {
private:
    TVector3 offset;
    TVector2 rotation;
    const Float_t Tel_DET_WIDTH;
    const Float_t Tel_DET_HEIGHT;
    const Float_t Tel_VSTRIP_WIDTH;
    const Float_t Tel_HSTRIP_WIDTH;
    
public:
    NewTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                        Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                        const TVector3& off, const TVector2& rot);
    
    bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      const std::vector<Float_t>& Edep_vec,
                      LightEjectile& ejectile) override;
};

// ========= PoP Detector Telescope Analyzer =========

class PoPTelescopeAnalyzer : public TelescopeAnalyzer {
private:
    const Float_t z_DE_pix;
    const Float_t Tel_DET_WIDTH;  // -10 to +10 mm
    const Float_t Tel_DET_HEIGHT;  // -10 to +10 mm
    const Float_t Tel_STRIP_WIDTH;
    const Float_t Tel_STRIP_OFFSET;
    const Float_t Tel_ROTATION_ANGLE;
    
    TTree* tree_E2;
    std::unordered_map<int, Float_t> E2_map;
    // Total energy reconstruction functions for each vertical strip
    // These functions reconstruct the true total energy from measured E_DE_E (E1 + DE)
    // Separate functions for low energy (not reaching E2) and high energy (reaching E2) events
    TF1* Etot_recon_strip_lowE[16];   // Low energy: events that stop in E1
    TF1* Etot_recon_strip_highE[16];  // High energy: events that reach E2
    
    void initializeEnergyReconstructionFunctions();
    
public:
    PoPTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                         Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                         TTree* te2);
    ~PoPTelescopeAnalyzer();
    
    bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      const std::vector<Float_t>& Edep_vec,
                      LightEjectile& ejectile) override;
};


// ========= Utility Functions =========

// Extract Z from PDGid
Int_t getZ(Int_t PDGid);

// Extract A from PDGid
Int_t getA(Int_t PDGid);

// Energy resolution function
Float_t getEnResolution(Float_t EnergyMeV, Float_t Res_Percent);

// GAGG scintillator resolution
double GAGG_resolution(double E);

// Read ejectile event file (Event_output text file) to get true values
std::unordered_map<int, std::tuple<double, double, double, int, int>> readEjectileFile(
    const char* reaction, const char* recType, Int_t excLabel, Double_t excEn,
    const ReactionInfo& reactionInfo);

// Helper function to load virtual detector data
struct VirtualDetectorData {
    std::unordered_map<int, std::tuple<Float_t, Float_t, Float_t>> pos_map;
    std::unordered_map<int, Int_t> pdgid_map;
};

VirtualDetectorData loadVirtualDetector(TFile* recoil_det_output, const char* tree_path,
                                        const char* reaction, const char* recType, 
                                        Int_t excLabel, Double_t excEn);

// Helper function to fill residue from virtual detector maps
void fillResidueFromMaps(Int_t eventID, 
                         const VirtualDetectorData& magsept_data,
                         const VirtualDetectorData& hrplane_data,
                         const VirtualDetectorData& quadwall_data,
                         HeavyResidue* residue);


// ========= Main analysis function =========
void det_analysis(Int_t excLabel, const char* recType);

#endif // DET_ANALYSIS_H
