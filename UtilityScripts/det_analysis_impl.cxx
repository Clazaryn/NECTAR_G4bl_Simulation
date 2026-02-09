// Implementation file for det_analysis.h
#include "UtilityScripts/det_analysis.h"

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

HeavyResidue::HeavyResidue() : Z(0), A(0), MagSept_x(0), MagSept_y(0), hit_MagSept(kFALSE),
                                HRplane_x(0), HRplane_y(0), hit_HRplane(kFALSE),
                                QuadWall_x(0), QuadWall_y(0), hit_QuadWall(kFALSE) {}

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

void PoPTelescopeAnalyzer::initializeEnergyReconstructionFunctions() {
    // Initialize total energy reconstruction functions for each vertical strip (1-16)
    // These functions reconstruct the true total energy Etot from measured E_DE_E (E1 + DE)
    // Separate functions for low energy (events stopping in E1) and high energy (events reaching E2)
    // Valid range: 3-18 MeV for E_DE_E
    // Coefficients from calibration fits to simulation data
    
    // Strip 1 (V01)
    Etot_recon_strip_lowE[0] = new TF1("Etot_recon_V01_lowE", "1.35069+0.978456*x-0.00489716*TMath::Power(x,2)+0.000226143*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[0] = new TF1("Etot_recon_V01_highE", "234.525-94.6461*x+19.1761*TMath::Power(x,2)-2.18803*TMath::Power(x,3)+0.143185*TMath::Power(x,4)-0.00500223*TMath::Power(x,5)+0.0000722282*TMath::Power(x,6)", 3, 18);
    
    // Strip 2 (V02)
    Etot_recon_strip_lowE[1] = new TF1("Etot_recon_V02_lowE", "1.3925+0.954073*x-0.00161447*TMath::Power(x,2)+0.000105024*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[1] = new TF1("Etot_recon_V02_highE", "233.788-94.1684*x+19.0509*TMath::Power(x,2)-2.17155*TMath::Power(x,3)+0.142*TMath::Power(x,4)-0.00496206*TMath::Power(x,5)+0.0000716648*TMath::Power(x,6)", 3, 18);
    
    // Strip 3 (V03)
    Etot_recon_strip_lowE[2] = new TF1("Etot_recon_V03_lowE", "1.20072+0.998161*x-0.00480596*TMath::Power(x,2)+0.000179519*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[2] = new TF1("Etot_recon_V03_highE", "231.546-92.742*x+18.6859*TMath::Power(x,2)-2.12404*TMath::Power(x,3)+0.138726*TMath::Power(x,4)-0.00484439*TMath::Power(x,5)+0.0000699967*TMath::Power(x,6)", 3, 18);
    
    // Strip 4 (V04)
    Etot_recon_strip_lowE[3] = new TF1("Etot_recon_V04_lowE", "1.21146+0.99131*x-0.00390203*TMath::Power(x,2)+0.000147161*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[3] = new TF1("Etot_recon_V04_highE", "230.306-92.0356*x+18.5171*TMath::Power(x,2)-2.1032*TMath::Power(x,3)+0.137345*TMath::Power(x,4)-0.00479821*TMath::Power(x,5)+0.0000693904*TMath::Power(x,6)", 3, 18);
    
    // Strip 5 (V05)
    Etot_recon_strip_lowE[4] = new TF1("Etot_recon_V05_lowE", "1.12839+1.01844*x-0.00644719*TMath::Power(x,2)+0.000221715*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[4] = new TF1("Etot_recon_V05_highE", "226.902-89.9009*x+17.978*TMath::Power(x,2)-2.03351*TMath::Power(x,3)+0.132476*TMath::Power(x,4)-0.00462336*TMath::Power(x,5)+0.0000668605*TMath::Power(x,6)", 3, 18);
    
    // Strip 6 (V06)
    Etot_recon_strip_lowE[5] = new TF1("Etot_recon_V06_lowE", "1.31035+0.965967*x-0.00128485*TMath::Power(x,2)+0.000055471*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[5] = new TF1("Etot_recon_V06_highE", "221.306-86.2871*x+17.0333*TMath::Power(x,2)-1.90566*TMath::Power(x,3)+0.12302*TMath::Power(x,4)-0.00426058*TMath::Power(x,5)+0.0000612131*TMath::Power(x,6)", 3, 18);
    
    // Strip 7 (V07)
    Etot_recon_strip_lowE[6] = new TF1("Etot_recon_V07_lowE", "1.40108+0.934843*x+0.00176821*TMath::Power(x,2)-0.0000340665*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[6] = new TF1("Etot_recon_V07_highE", "211.533-79.5878*x+15.1947*TMath::Power(x,2)-1.64686*TMath::Power(x,3)+0.10327*TMath::Power(x,4)-0.0034838*TMath::Power(x,5)+0.0000488816*TMath::Power(x,6)", 3, 18);
    
    // Strip 8 (V08)
    Etot_recon_strip_lowE[7] = new TF1("Etot_recon_V08_lowE", "1.49675+0.899169*x+0.005326*TMath::Power(x,2)-0.00014451*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[7] = new TF1("Etot_recon_V08_highE", "181.925-58.2245*x+9.0296*TMath::Power(x,2)-0.736542*TMath::Power(x,3)+0.0306613*TMath::Power(x,4)-0.000511943*TMath::Power(x,5)", 3, 18);
    
    // Strip 9 (V09)
    Etot_recon_strip_lowE[8] = new TF1("Etot_recon_V09_lowE", "1.64485+0.849583*x+0.0101998*TMath::Power(x,2)-0.000290011*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[8] = new TF1("Etot_recon_V09_highE", "179.959025-57.2111*x+8.829424*TMath::Power(x,2)-0.717739*TMath::Power(x,3)+0.029820676*TMath::Power(x,4)-0.000497579*TMath::Power(x,5)", 3, 18);
    
    // Strip 10 (V10)
    Etot_recon_strip_lowE[9] = new TF1("Etot_recon_V10_lowE", "1.69229+0.842741*x+0.0106049*TMath::Power(x,2)-0.000302836*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[9] = new TF1("Etot_recon_V10_highE", "177.335-55.7246*x+8.50325*TMath::Power(x,2)-0.683296*TMath::Power(x,3)+0.028073*TMath::Power(x,4)-0.000463473*TMath::Power(x,5)", 3, 18);
    
    // Strip 11 (V11)
    Etot_recon_strip_lowE[10] = new TF1("Etot_recon_V11_lowE", "1.70029+0.843103*x+0.0106755*TMath::Power(x,2)-0.000313177*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[10] = new TF1("Etot_recon_V11_highE", "174.317-54.0509*x+8.15242*TMath::Power(x,2)-0.648289*TMath::Power(x,3)+0.026403*TMath::Power(x,4)-0.000432896*TMath::Power(x,5)", 3, 18);
    
    // Strip 12 (V12)
    Etot_recon_strip_lowE[11] = new TF1("Etot_recon_V12_lowE", "1.76783+0.821267*x+0.0126448*TMath::Power(x,2)-0.000364064*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[11] = new TF1("Etot_recon_V12_highE", "170.735-52.0978*x+7.748*TMath::Power(x,2)-0.608218*TMath::Power(x,3)+0.0244946*TMath::Power(x,4)-0.000397816*TMath::Power(x,5)", 3, 18);
    
    // Strip 13 (V13)
    Etot_recon_strip_lowE[12] = new TF1("Etot_recon_V13_lowE", "1.92214+0.76626*x+0.0182469*TMath::Power(x,2)-0.000536454*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[12] = new TF1("Etot_recon_V13_highE", "167.925-50.6538*x+7.46225*TMath::Power(x,2)-0.58078*TMath::Power(x,3)+0.023211*TMath::Power(x,4)-0.000374337*TMath::Power(x,5)", 3, 18);
    
    // Strip 14 (V14)
    Etot_recon_strip_lowE[13] = new TF1("Etot_recon_V14_lowE", "1.93006+0.758802*x+0.0193614*TMath::Power(x,2)-0.000581656*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[13] = new TF1("Etot_recon_V14_highE", "169.777-51.6241*x+7.65851*TMath::Power(x,2)-0.599727*TMath::Power(x,3)+0.0240867*TMath::Power(x,4)-0.000389921*TMath::Power(x,5)", 3, 18);
    
    // Strip 15 (V15)
    Etot_recon_strip_lowE[14] = new TF1("Etot_recon_V15_lowE", "1.94292+0.754089*x+0.0197303*TMath::Power(x,2)-0.00058678*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[14] = new TF1("Etot_recon_V15_highE", "163.793-48.3521*x+6.97445*TMath::Power(x,2)-0.531018*TMath::Power(x,3)+0.020763*TMath::Power(x,4)-0.000327826*TMath::Power(x,5)", 3, 18);
    
    // Strip 16 (V16)
    Etot_recon_strip_lowE[15] = new TF1("Etot_recon_V16_lowE", "1.93858+0.75949*x+0.0188402*TMath::Power(x,2)-0.000548706*TMath::Power(x,3)", 3, 18);
    Etot_recon_strip_highE[15] = new TF1("Etot_recon_V16_highE", "158.935-45.8074*x+6.46525*TMath::Power(x,2)-0.482008*TMath::Power(x,3)+0.0184859*TMath::Power(x,4)-0.000286842*TMath::Power(x,5)", 3, 18);
}

PoPTelescopeAnalyzer::PoPTelescopeAnalyzer(const char* r, const char* rt, Int_t el,
                         Double_t mb, Double_t mt, Double_t mr, Double_t me, Double_t ek, Double_t p,
                         TTree* te2)
    : TelescopeAnalyzer(r, rt, el, mb, mt, mr, me, ek, p), tree_E2(te2) {
    
    initializeEnergyReconstructionFunctions();
    
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
        delete Etot_recon_strip_lowE[i];
        delete Etot_recon_strip_highE[i];
    }
}

bool PoPTelescopeAnalyzer::analyzeEvent(Int_t eventID, Float_t x_DE, Float_t y_DE, Float_t z_DE,
                      Float_t Edep_DE, Float_t Edep_E1, Float_t Edep_Eres,
                      LightEjectile& ejectile) {
    
    // Calculate strip numbers (PoP geometry)
    // Strip numbering is inverted: strip 1 is at x=10 (rightmost), strip 16 is at x=-10 (leftmost)
    // Original formula: Tel_VSTRIP = 17 - ceil((x_DE + 10) / 1.25)
    Int_t vert_strip = 17 - ceil((x_DE + Tel_STRIP_OFFSET) / Tel_STRIP_WIDTH);
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
    // Apply detector resolution to measured energies
    Float_t DSSSD = getEres(Edep_DE, 1.0);
    Float_t Etot_E = getEres(Edep_E1, 1.1);
    Float_t E_DE_E = Etot_E + DSSSD;  // Measured energy: E1 + DE
    
    // Reconstruct true total energy using strip-specific and energy-dependent functions
    // Low energy: events that stop in E1 (not reaching E2)
    // High energy: events that reach E2 (penetrating through E1)
    Float_t Etot_tel_res;
    if (isInE2) {
        Etot_tel_res = Etot_recon_strip_highE[vert_strip-1]->Eval(E_DE_E);
    } else {
        Etot_tel_res = Etot_recon_strip_lowE[vert_strip-1]->Eval(E_DE_E);
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

// ========= Helper Functions Implementation =========

VirtualDetectorData loadVirtualDetector(TFile* recoil_file, const char* tree_path) {
    VirtualDetectorData data;
    
    TTree* tree = (TTree*)recoil_file->Get(tree_path);
    if (!tree) {
        return data;  // Return empty data if tree not found
    }
    
    Float_t EventID, x, y, z, PDGid;
    tree->SetBranchAddress("EventID", &EventID);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    
    bool has_pdgid = (tree->GetBranch("PDGid") != nullptr);
    if (has_pdgid) {
        tree->SetBranchAddress("PDGid", &PDGid);
    }
    
    for (Int_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        Int_t eventID = static_cast<int>(EventID);
        data.pos_map[eventID] = std::make_tuple(x, y, z);
        if (has_pdgid) {
            data.pdgid_map[eventID] = PDGid;
        }
    }
    
    return data;
}

void fillResidueFromMaps(Int_t eventID, 
                         const VirtualDetectorData& magsept_data,
                         const VirtualDetectorData& hrplane_data,
                         const VirtualDetectorData& quadwall_data,
                         HeavyResidue* residue) {
    // Fill MagSept data if available
    if (magsept_data.pos_map.count(eventID) > 0) {
        auto mag_pos = magsept_data.pos_map.at(eventID);
        residue->MagSept_x = std::get<0>(mag_pos);
        residue->MagSept_y = std::get<1>(mag_pos);
        residue->hit_MagSept = kTRUE;
        
        // Extract Z and A from PDGid if available
        if (magsept_data.pdgid_map.count(eventID) > 0) {
            Float_t pdgid = magsept_data.pdgid_map.at(eventID);
            residue->Z = getZ(pdgid);
            residue->A = getA(pdgid);
        }
    } else {
        residue->hit_MagSept = kFALSE;
    }
    
    // Fill HRplane data if available
    if (hrplane_data.pos_map.count(eventID) > 0) {
        auto hr_pos = hrplane_data.pos_map.at(eventID);
        residue->HRplane_x = std::get<0>(hr_pos);
        residue->HRplane_y = std::get<1>(hr_pos);
        residue->hit_HRplane = kTRUE;
    } else {
        // Event hit wall before reaching HRplane
        residue->HRplane_x = 0;
        residue->HRplane_y = 0;
        residue->hit_HRplane = kFALSE;
    }
    
    // Fill QuadWall data if available
    if (quadwall_data.pos_map.count(eventID) > 0) {
        auto qw_pos = quadwall_data.pos_map.at(eventID);
        residue->QuadWall_x = std::get<0>(qw_pos);
        residue->QuadWall_y = std::get<1>(qw_pos);
        residue->hit_QuadWall = kTRUE;
    } else {
        // Event did not hit QuadWall
        residue->QuadWall_x = 0;
        residue->QuadWall_y = 0;
        residue->hit_QuadWall = kFALSE;
    }
}
