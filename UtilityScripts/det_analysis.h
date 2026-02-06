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

// ========= Utility Functions =========

// Extract Z from PDGid
// PDGid format: "100" + Z_out (3 digits) + A_out (3 digits) + "0"
// Example: Z=1, A=1 -> "1000010010" -> 1000010010
// Z extraction: (PDGid - 1000000000) / 10000, then floor
Int_t getZ(Float_t PDGid);

// Extract A from PDGid
// A extraction: (PDGid - 1000000000 - Z*10000) / 10
Int_t getA(Float_t PDGid);

// Energy resolution function
Float_t getEres(Float_t EnergyMeV, Float_t Res_Percent);

// GAGG scintillator resolution
double GAGG_resolution(double E);

// Read ejectile file to get true values
std::unordered_map<int, std::tuple<double, double, double, int, int>> readEjectileFile(
    const char* reaction, const char* recType, Int_t excLabel);

// ========= Data Classes =========

class LightEjectile {
public:
    Int_t Z, A;
    Double_t true_Estar, true_Etot, true_theta;
    Double_t recon_Estar, recon_Etot, recon_theta;
    Int_t vert_strip;
    Int_t detector_id;  // 0=primary, 1=auxillary for New detectors
    
    LightEjectile();
};

class HeavyResidue {
public:
    Int_t Z, A;
    Double_t MagSept_x, MagSept_y;
    Bool_t has_MagSept;  // true if event has MagSept data
    Double_t HRplane_x, HRplane_y;
    Bool_t has_HRplane;  // true if event has HRplane data
    
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
                              Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                              LightEjectile& ejectile) = 0;
};

// ========= New Detector Telescope Analyzer =========

class NewTelescopeAnalyzer : public TelescopeAnalyzer {
private:
    TVector3 offset;
    TVector2 rotation;
    static const Float_t Tel_DET_WIDTH;
    static const Float_t Tel_DET_HEIGHT;
    static const Float_t Tel_VSTRIP_WIDTH;
    static const Float_t Tel_HSTRIP_WIDTH;
    
public:
    NewTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                        Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                        const TVector3& off, const TVector2& rot);
    
    bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) override;
};

// ========= PoP Detector Telescope Analyzer =========

class PoPTelescopeAnalyzer : public TelescopeAnalyzer {
private:
    static const Float_t z_DE_pix;
    static const Float_t Tel_DET_WIDTH;
    static const Float_t Tel_DET_HEIGHT;
    static const Float_t Tel_STRIP_WIDTH;
    static const Float_t Tel_STRIP_OFFSET;
    static const Float_t Tel_ROTATION_ANGLE;
    
    TTree* tree_E2;
    std::unordered_map<int, Float_t> E2_map;
    TF1* E_corr_lowE_prot[16];
    TF1* E_corr_highE_prot[16];
    
    void initializeCorrectionFunctions();
    
public:
    PoPTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                         Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                         TTree* te2);
    ~PoPTelescopeAnalyzer();
    
    bool analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) override;
};

#endif // DET_ANALYSIS_H
