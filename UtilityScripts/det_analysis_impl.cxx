// Implementation file for det_analysis.h
#include "CppHeaders/det_analysis.h"

// ========= Utility Functions Implementation =========

Int_t getZ(Float_t PDGid) {
    // PDGid format: "100" + Z_out (3 digits) + A_out (3 digits) + "0"
    // Extract Z: (PDGid - 1000000000) / 10000, then floor
    return floor((PDGid - 1000000000) * 0.0001);
}

Int_t getA(Float_t PDGid) {
    // Extract A: (PDGid - 1000000000 - Z*10000) / 10
    Int_t Z = getZ(PDGid);
    return (PDGid - 1000000000 - Z * 10000) * 0.1;
}

Float_t getEres(Float_t EnergyMeV, Float_t Res_Percent) {
    return gRandom->Gaus(EnergyMeV, Res_Percent / 100 * EnergyMeV);
}

double GAGG_resolution(double E) {
    const double a = 11.2262, b = -0.573112;
    if (E <= 0) return -1.0;
    else if (E <= 1) return 10.0 / 2.355;
    else return (a * pow(E, b)) / 2.355;
}

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

// ========= Data Classes Implementation =========

LightEjectile::LightEjectile() : Z(0), A(0), true_Estar(0), true_Etot(0), true_theta(0),
                                  recon_Estar(0), recon_Etot(0), recon_theta(0), vert_strip(0), detector_id(0) {}

HeavyResidue::HeavyResidue() : Z(0), A(0), MagSept_x(0), MagSept_y(0), has_MagSept(kFALSE),
                                HRplane_x(0), HRplane_y(0), has_HRplane(kFALSE) {}

// ========= Telescope Analyzer Implementation =========

TelescopeAnalyzer::TelescopeAnalyzer(const char* r, const char* rt, Int_t el, 
                                     Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p) 
    : reaction(r), recType(rt), excLabel(el), mass_beam(mb), mass_target(mt), 
      mass_recoil(mr), mass_ejectile(me), Ek_beam(ek), P_beam(p) {}

// ========= New Detector Telescope Analyzer Implementation =========

const Float_t NewTelescopeAnalyzer::Tel_DET_WIDTH = 122.0;
const Float_t NewTelescopeAnalyzer::Tel_DET_HEIGHT = 40.0;
const Float_t NewTelescopeAnalyzer::Tel_VSTRIP_WIDTH = 1.0;
const Float_t NewTelescopeAnalyzer::Tel_HSTRIP_WIDTH = 1.0;

NewTelescopeAnalyzer::NewTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                        Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                        const TVector3& off, const TVector2& rot)
    : TelescopeAnalyzer(r, rt, el, mb, mt, mr, me, ek, p), offset(off), rotation(rot) {}

bool NewTelescopeAnalyzer::analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) {
    
    // Calculate vertical strip number
    Int_t vert_strip = ceil((x_DE + Tel_DET_WIDTH/2) / Tel_VSTRIP_WIDTH);
    Int_t hor_strip = ceil((y_DE + Tel_DET_HEIGHT/2) / Tel_HSTRIP_WIDTH);
    
    // throw away events in the edge strips of the DSSD to avoid edge effects
    if (!(1 < hor_strip && hor_strip < Tel_DET_HEIGHT/Tel_HSTRIP_WIDTH) || !(1 < vert_strip && vert_strip < Tel_DET_WIDTH/Tel_VSTRIP_WIDTH)) {
        return false;
    }
    
    // Calculate pixel center coordinates
    Float_t x_DE_pix = -Tel_DET_WIDTH/2 + Tel_VSTRIP_WIDTH*(vert_strip-1) + Tel_VSTRIP_WIDTH/2;
    Float_t y_DE_pix = -Tel_DET_HEIGHT/2 + Tel_HSTRIP_WIDTH*(hor_strip-1) + Tel_HSTRIP_WIDTH/2;
    
    // Apply rotation to get absolute coordinates
    // Intrinsic rotations: rotate by alpha about y and then by beta about x
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

// ========= PoP Detector Telescope Analyzer Implementation =========

const Float_t PoPTelescopeAnalyzer::z_DE_pix = 99.025;
const Float_t PoPTelescopeAnalyzer::Tel_DET_WIDTH = 20.0;  // -10 to +10 mm
const Float_t PoPTelescopeAnalyzer::Tel_DET_HEIGHT = 20.0;  // -10 to +10 mm
const Float_t PoPTelescopeAnalyzer::Tel_STRIP_WIDTH = 1.25;
const Float_t PoPTelescopeAnalyzer::Tel_STRIP_OFFSET = 10.0;
const Float_t PoPTelescopeAnalyzer::Tel_ROTATION_ANGLE = 60.0 * M_PI / 180.0;

void PoPTelescopeAnalyzer::initializeCorrectionFunctions() {
    // Initialize energy correction functions for each vertical strip
    // These are simplified - full implementation would include all 16 strips with exact coefficients
    // For now, using placeholder functions
    for (int i = 0; i < 16; ++i) {
        E_corr_lowE_prot[i] = new TF1(Form("E_corr_V%02d_lowE_prot", i+1),
            "1.3+0.98*x-0.005*TMath::Power(x,2)+0.0002*TMath::Power(x,3)", 3, 18);
        E_corr_highE_prot[i] = new TF1(Form("E_corr_V%02d_highE_prot", i+1),
            "200-80*x+16*TMath::Power(x,2)-2*TMath::Power(x,3)+0.14*TMath::Power(x,4)-0.005*TMath::Power(x,5)", 3, 18);
    }
}

PoPTelescopeAnalyzer::PoPTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
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

PoPTelescopeAnalyzer::~PoPTelescopeAnalyzer() {
    for (int i = 0; i < 16; ++i) {
        delete E_corr_lowE_prot[i];
        delete E_corr_highE_prot[i];
    }
}

bool PoPTelescopeAnalyzer::analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) {
    
    // Calculate strip numbers (PoP geometry)
    // Strip calculation: Tel_VSTRIP = 17 - ceil((x_DE + 10) / 1.25)
    Int_t vert_strip = ceil((x_DE + Tel_STRIP_OFFSET) / Tel_STRIP_WIDTH);
    Int_t hor_strip = ceil((y_DE + Tel_STRIP_OFFSET) / Tel_STRIP_WIDTH);
    
    // throw away events in the edge strips of the DSSD to avoid edge effects
    if (!(1 < hor_strip && hor_strip < Tel_DET_HEIGHT/Tel_STRIP_WIDTH) || !(1 < vert_strip && vert_strip < Tel_DET_WIDTH/Tel_STRIP_WIDTH)) {
        return false;
    }
    
    // Calculate pixel center coordinates
    Float_t x_DE_pix = -Tel_STRIP_OFFSET + Tel_STRIP_WIDTH*(vert_strip-1) + Tel_STRIP_WIDTH/2;
    Float_t y_DE_pix = -Tel_STRIP_OFFSET + Tel_STRIP_WIDTH*(hor_strip-1) + Tel_STRIP_WIDTH/2;
    
    // PoP theta calculation with rotation
    // Apply 60 degree rotation: zz = z*cos(60°) + x*cos(30°), xx = z*sin(60°) - x*sin(30°)
    Float_t zz = z_DE_pix*cos(Tel_ROTATION_ANGLE) + x_DE_pix*cos(M_PI/6);  // cos(30°) = cos(π/6)
    Float_t xx = z_DE_pix*sin(Tel_ROTATION_ANGLE) - x_DE_pix*sin(M_PI/6);   // sin(30°) = sin(π/6)
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
        Etot_tel_res = E_corr_highE_prot[vert_strip-1]->Eval(E_DE_E);
    } else {
        Etot_tel_res = E_corr_lowE_prot[vert_strip-1]->Eval(E_DE_E);
    }
    
    // Calculate excitation energy
    Double_t P_tel_res = sqrt(pow(Etot_tel_res + mass_ejectile, 2) - pow(mass_ejectile, 2));
    Double_t Exc_energy = sqrt(pow(Ek_beam + mass_beam + mass_target - Etot_tel_res - mass_ejectile, 2.) -
                               pow(P_beam, 2.) - pow(P_tel_res, 2.) + 2*P_beam*P_tel_res*TMath::Cos(theta_DE_pix)) - mass_recoil;
    
    // Fill ejectile structure
    ejectile.vert_strip = vert_strip;
    ejectile.recon_theta = theta_DE_pix * 180.0 / M_PI;
    ejectile.recon_Etot = Etot_tel_res;
    ejectile.recon_Estar = Exc_energy;
    
    return true;
}
