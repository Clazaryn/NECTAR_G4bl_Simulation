// --
// NECTAR simulation of new chamber for 206Pb experiment
// Last modification : G. Leckenby - 1/07/2025
// Adapted from ExcResolution code from C. Berthelot

//  ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
//
//      	Simulation analysis 
//
//  ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

// ========= C++ libraries =========
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <unordered_map>
#include <sstream>
using namespace std;

// ========= ROOT libraries =========
#include <TROOT.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TGraph.h>
#include <TSpline.h>

// ========= AI header files =========
#include "UtilityScripts/reaction_info.h"		// Include reaction information parser

// ========= Variable declaration ==========
// Global variables for reaction information
ReactionInfo reactionInfo;
const char* reaction;
Double_t mass_beam, mass_target, mass_recoil, mass_ejectile;
Float_t EuA_beam;
Double_t Ek_beam, P_beam;

// Common parameters
const Float_t   Tel_DET_WIDTH=122.,  Tel_DET_HEIGHT=40.;
const Float_t 	Tel_VSTRIP_WIDTH=0.01, Tel_HSTRIP_WIDTH=0.01;

// Front pocket parameters
const Float_t 	FPT_X_Rotation=-20*M_PI/180, FPT_Y_Rotation=-20*M_PI/180;
const Float_t 	FPT_X_Offset=(-80.858+3.025*sin(FPT_Y_Rotation)*cos(FPT_X_Rotation))+0*sin(FPT_Y_Rotation);
const Float_t 	FPT_Y_Offset=(49.802+3.025*sin(FPT_X_Rotation));
const Float_t 	FPT_Z_Offset=(31.162+3.025*cos(FPT_Y_Rotation)*cos(FPT_X_Rotation))+0*cos(FPT_Y_Rotation);

// Back pocket parameters
const Float_t 	BPT_X_Rotation=-20*M_PI/180, BPT_Y_Rotation=0;
const Float_t 	BPT_X_Offset=(0+3.025*sin(BPT_Y_Rotation)*cos(BPT_X_Rotation));
const Float_t 	BPT_Y_Offset=(49.802+3.025*sin(BPT_X_Rotation));
const Float_t 	BPT_Z_Offset=(206.54+3.025*cos(BPT_Y_Rotation)*cos(BPT_X_Rotation));

const std::vector<TVector3> offsets = {
	TVector3(FPT_X_Offset, FPT_Y_Offset, FPT_Z_Offset),
	TVector3(BPT_X_Offset, BPT_Y_Offset, BPT_Z_Offset)
};	
const std::vector<TVector2> rotations = {
	TVector2(FPT_X_Rotation, FPT_Y_Rotation),
	TVector2(BPT_X_Rotation, BPT_Y_Rotation)
};

// ========= Functions =========
// Initialize reaction information from INI file
void initializeReactionInfo() {
    if (!reactionInfo.loadFromIni("reac_info_nChamb.txt")) {
        std::cerr << "Error: Could not load reaction info from INI file. Exiting." << std::endl;
        return;
    }
    
    // Set global variables from INI file
    reaction = reactionInfo.reaction.c_str();
    mass_beam = reactionInfo.mass_beam;
    mass_target = reactionInfo.mass_target;
    mass_recoil = reactionInfo.mass_recoil;
    mass_ejectile = reactionInfo.mass_eject;
    EuA_beam = reactionInfo.beam_EuA;
    Ek_beam = mass_beam / 931.49410242 * EuA_beam;	 
    P_beam = sqrt(pow(Ek_beam + mass_beam, 2.0) - pow(mass_beam, 2.0));
}

// Adds telescope energy resolution gaussian to energy values
Float_t getEres(Float_t EnergyMeV, Float_t Res_Percent)
{	
	Float_t E_res = gRandom->Gaus(EnergyMeV,Res_Percent/100*EnergyMeV);	
	return E_res;
}

// determines if file exists
bool file_exists(const char* filename) 
{
    std::ifstream infile(filename);
    return infile.good();
}

// function to open or create root file
TFile* open_or_create_rootfile(const char* filename) {
    if (file_exists(filename)) {
        std::cout << "Opening existing file: " << filename << std::endl;
        TFile* f = TFile::Open(filename, "UPDATE");
        if (!f || f->IsZombie()) { std::cerr << "Error opening file: " << filename << std::endl; return nullptr; }
        return f;
    } else {
        std::cout << "Creating new file: " << filename << std::endl;
        TFile* f = new TFile(filename, "RECREATE");
        if (!f || f->IsZombie()) { std::cerr << "Error creating file: " << filename << std::endl; return nullptr; }
        return f;
    }
}

// creates histogram labels from floats
TString labelFromFloat(const TString& prefix, Float_t x, const TString& suffix = "", Int_t precision = 2) 
{
    char floatStr[50];
    sprintf(floatStr, "%.*f", precision, x);  // format float with precision

    // Replace '.' with 'p'
    char* dot = strchr(floatStr, '.');
    if (dot) *dot = 'p';

    // Construct final label
    TString label = prefix + floatStr + suffix;
    return label;
}

// power law description of GAGG scintillator resolution
double GAGG_resolution(double E) 
{
	// constants come from power law fit to two points from Furano et al (2021) JINST 16:P10012
    const double a = 11.2262, b = -0.573112;
    if (E <= 0) { return -1.0; }				// avoid invalid power domain
	else if (E<= 1)	{ return 10/2.355; }		// minimum 10% resolution to avoid E->0 blowup
	else { return (a * pow(E, b))/2.355; }		// returns in sigma but calibrated on FWHM
}

// interpolating function of energy lost in SS window for proton of varying energy
double pInSS_interp(double x_query) 
{
    // Below data is a result of SRIM simulation of proton in Stainless steel window
	// x_vals are energy after window (MeV) and y_vals are energy deposited in window (MeV)
    const int n = 11;
    double x_vals[n] = {0.32084, 1.36165, 2.77431, 3.9905, 6.77083, 9.41821, 14.5745, 24.7129, 49.8315, 99.899, 149.924};
    double y_vals[n] = {2.17916, 1.63835, 1.22569, 1.0095, 0.729167, 0.581787, 0.425458, 0.2871, 0.168498, 0.100955, 0.0755165};

    TGraph* graph = new TGraph(n, x_vals, y_vals);			// Create TGraph from the data points
	TSpline3* spline = new TSpline3("spline", graph);		// Create a cubic spline from the graph
	double result = spline->Eval(x_query);					// Evaluate the spline at the desired x value

    delete spline;    delete graph;
    return result;
}

// function to read ejectile file and create map of EventID to true values
std::unordered_map<int, std::pair<double, double>> readEjectileFile(const char* reaction, const char* recType, Int_t excLabel) {
    std::unordered_map<int, std::pair<double, double>> eventMap;
    
    TString filename = Form("../%s_sim/Event_output/output_event_generator_%s_%s_excEn%02d_ejectile.txt", 
                           reaction, reaction, recType, excLabel);
    
    std::ifstream file(filename.Data());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open ejectile file: " << filename << std::endl;
        return eventMap;
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x, y, z, px, py, pz, t;
        int pdgid, eventid, trackid, parentid;
        double weight;
        
        if (iss >> x >> y >> z >> px >> py >> pz >> t >> pdgid >> eventid >> trackid >> parentid >> weight) {
            // Calculate true theta from momentum components
            double p_mag = sqrt(px*px + py*py + pz*pz);
            double theta_true = acos(pz/p_mag) * 180.0 / M_PI;
            
            // Calculate true kinetic energy (assuming proton mass)
            double proton_mass = 938.7830736444523; // MeV/c^2
            double E_true = sqrt(p_mag*p_mag + proton_mass*proton_mass) - proton_mass;
            
            eventMap[eventid] = std::make_pair(theta_true, E_true);
        }
    }
    
    file.close();
    std::cout << "Loaded " << eventMap.size() << " events from ejectile file" << std::endl;
    return eventMap;
}


// function to analyse single telescope detector and fill histograms
void analyseTelescope( const char* reaction, const char* recType, Int_t excLabel, const TString& name,
    TTree* tree_DE, TTree* tree_E1, TTree* tree_Eres, std::vector<TTree*> trees_E1_veto, std::vector<TTree*> trees_Eres_veto, TTree* Dip_block, TTree* HR_block,
    const TVector3& offset, const TVector2& rotation, const Int_t nThetaBins, const Double_t thetaBinEdges[],
    TH1D* hEex, TH1D* hTheta, TH1D* hEprot, TH1D* hThetaAccur, TH1D* hEprotAccur, TH2D* hThetaAccurvsTheta, TH2D* hEprotAccurvsEprot,
	TH2D* hDEvE1Ban, TH2D* hDEvEresBan, TH2D* hE1Veto_dE_posMap, TH2D* hEresVeto_dE_posMap, std::vector<TH2D*> hDEvE1Ban_thetaBins, std::vector<TH2D*> hDEvEresBan_thetaBins,
	TH2D* hDipCoinc, TH2D* hHRcoinc, TH2D* hHRvEexcCoinc, std::vector<TH2D*> hHRvEexcCoinc_thetaBins
) {
   	std::unordered_map<int, std::pair<double, double>> eventMap = readEjectileFile(reaction, recType, excLabel);	// Read ejectile file to get true values for accuracy measurement
   	
   	// ========= Variables for analysis ==========
	Float_t En_tot, theta_calc, P_ejc, Exc_en;
   	// ========= Variables for tree branches ==========
   	Float_t nEntries_DE, x_DE, y_DE, z_DE, EventID_DE, Edep_DE, Erec_DE=0, Erec_DE_bl=0;		// Telescope - DeltaE detector
   	Float_t nEntries_E1, EventID_E1, Edep_E1, Erec_E1=0, Erec_E1_bl=0;							// Telescope - E1 detector
	Float_t nEntries_Eres, EventID_Eres, Edep_Eres, Erec_Eres=0, Erec_Eres_bl=0;				// Telescope - Eres detector
	Float_t nEntries_Dip, x_Dip, y_Dip, z_Dip, EventID_Dip, Edep_Dip, Erec_Dip=0;				// Dipole entrance detector
	Float_t nEntries_HR, x_HR, y_HR, z_HR, EventID_HR, Edep_HR, Erec_HR=0;						// Heavy Residue detector

	
   	// ========= Tree branching ==========
	// Telescope - DeltaE detector
  	nEntries_DE = tree_DE->GetEntries(); 				
  	tree_DE->SetBranchAddress("EventID", &EventID_DE);		tree_DE->SetBranchAddress("x", &x_DE);			
	tree_DE->SetBranchAddress("y", &y_DE);					tree_DE->SetBranchAddress("z", &z_DE);		
	tree_DE->SetBranchAddress("Edep", &Edep_DE);
	// Telescope - E1 detector
	nEntries_E1 = tree_E1->GetEntries();
	tree_E1->SetBranchAddress("Edep", &Edep_E1);			tree_E1->SetBranchAddress("EventID", &EventID_E1);
	// Telescope - Eres detector
	nEntries_Eres = tree_Eres->GetEntries();			
   	tree_Eres->SetBranchAddress("Edep", &Edep_Eres);   		tree_Eres->SetBranchAddress("EventID", &EventID_Eres);
	// Telescope - Dipole entrance block detector
	nEntries_Dip = Dip_block->GetEntries(); 				
	Dip_block->SetBranchAddress("EventID", &EventID_Dip);		Dip_block->SetBranchAddress("x", &x_Dip);			
  	Dip_block->SetBranchAddress("y", &y_Dip);					Dip_block->SetBranchAddress("z", &z_Dip);
	// Telescope - HR plane block detector
	nEntries_HR = HR_block->GetEntries(); 				
	HR_block->SetBranchAddress("EventID", &EventID_HR);		HR_block->SetBranchAddress("x", &x_HR);			
  	HR_block->SetBranchAddress("y", &y_HR);					HR_block->SetBranchAddress("z", &z_HR);

	// Telescope - veto trees
	std::vector<Float_t> nEntries_E1Veto(4), EventID_E1Veto(4), nEntries_EresVeto(4), EventID_EresVeto(4);
	for (int i = 0; i < 4; ++i) {
		nEntries_E1Veto[i] = trees_E1_veto[i]->GetEntries();		trees_E1_veto[i]->SetBranchAddress("EventID", &EventID_E1Veto[i]);
		nEntries_EresVeto[i] = trees_Eres_veto[i]->GetEntries();	trees_Eres_veto[i]->SetBranchAddress("EventID", &EventID_EresVeto[i]);
	}

	// ========= Build lookup maps for EventID matching ========= - avoids recursive for loops
	std::unordered_map<int, Float_t> E1_map, Eres_map;
	std::vector<std::unordered_set<int>> E1Veto_list(4), EresVeto_list(4);
	std::unordered_map<int, std::tuple<Float_t,Float_t,Float_t>> Dip_map, HR_map;

	for (int i = 0; i < nEntries_E1; ++i) {
	    tree_E1->GetEntry(i);
	    E1_map[static_cast<int>(EventID_E1)] = Edep_E1;  // Only keeps the last if duplicates
	}

	for (int i = 0; i < nEntries_Eres; ++i) {
	    tree_Eres->GetEntry(i);
	    Eres_map[static_cast<int>(EventID_Eres)] = Edep_Eres;
	}

	for (int j = 0; j < 4; ++j) {
		for (int i = 0; i < nEntries_E1Veto[j]; ++i) {
			trees_E1_veto[j]->GetEntry(i);		E1Veto_list[j].insert(static_cast<int>(EventID_E1Veto[j]));
		}
		for (int i = 0; i < nEntries_EresVeto[j]; ++i) {
			trees_Eres_veto[j]->GetEntry(i);	EresVeto_list[j].insert(static_cast<int>(EventID_EresVeto[j]));
		}
	}
	
	for (int i = 0; i < nEntries_Dip; ++i) {
	    Dip_block->GetEntry(i);
	    Dip_map[static_cast<int>(EventID_Dip)] = std::make_tuple(x_Dip, y_Dip, z_Dip);
	}

	for (int i = 0; i < nEntries_HR; ++i) {
	    HR_block->GetEntry(i);
	    HR_map[static_cast<int>(EventID_HR)] = std::make_tuple(x_HR, y_HR, z_HR);
	}
   	
	// ========= Event loop gating on DSSD ==========
  	for(Int_t nevent_DE=0; nevent_DE<nEntries_DE; nevent_DE++) {
  		tree_DE->GetEntry(nevent_DE);
		int EventID_cast = static_cast<int>(EventID_DE);		// static cast required because floats could mess with map lookup
		//std::cout << "EventID_cast is " << EventID_cast << std::endl;

  		//Telescope strip calculation from ejectile coordinates
		// x/y dimensions are given wrt centre of 122x40 DSSD - assuming 1 mm strips
  		Int_t vert_strip = ceil((x_DE+Tel_DET_WIDTH/2)/Tel_VSTRIP_WIDTH);		
		Int_t hor_strip	 = ceil((y_DE+Tel_DET_HEIGHT/2)/Tel_HSTRIP_WIDTH);
		//cout << "VSTRIP is: " << Tel_VSTRIP << ", HSTRIP is: " << Tel_HSTRIP << endl;

  		if (1<hor_strip && hor_strip<40/Tel_HSTRIP_WIDTH) {
	  		// Telescope angle calculation pixel coordinate (taken at center of each pixel)
			// need to extract absolute x,y,z coordinates of pixel centre. x/y_DE_pix stores relative position
			Float_t x_DE_pix = -Tel_DET_WIDTH/2 + Tel_VSTRIP_WIDTH*(vert_strip-1) + Tel_VSTRIP_WIDTH/2;	
			Float_t y_DE_pix = -Tel_DET_HEIGHT/2 + Tel_HSTRIP_WIDTH*(hor_strip-1) + Tel_HSTRIP_WIDTH/2;
			// to determine absolute, we need to account for rotation of new telescope about two axes
			// we apply intrinsic rotations: rotate by alpha about y and then by beta about x. 
			// intrinsic rotation means rotation matrix is Ry(alpha).Rx(beta).v_DE
			// Ry(alpha) = {{Cos[alpha], 0, Sin[alpha]}, {0, 1, 0}, {-Sin[alpha], 0, Cos[alpha]}}
			// Rx(beta) = {{1, 0, 0}, {0, Cos[beta], -Sin[beta]}, {0, Sin[beta], Cos[beta]}}
			Float_t x_DE_abs = offset.X() + (x_DE_pix*cos(rotation.Y()) + y_DE_pix*sin(rotation.Y())*sin(rotation.X()));
			Float_t y_DE_abs = offset.Y() + (y_DE_pix*cos(rotation.X()));
			Float_t z_DE_abs = offset.Z() + (-x_DE_pix*sin(rotation.Y()) + y_DE_pix*cos(rotation.Y())*sin(rotation.X()));

			// it is now straightforward to calculate theta from standard spherical coordinates
			TVector3 abs_vec(x_DE_abs, y_DE_abs, z_DE_abs);
			theta_calc = abs_vec.Theta() * 180.0 / M_PI;

			// assign energy values for detection
			Erec_DE = Edep_DE;
			// Look up matching EventID from E1 and Eres - this is a ternary query based on our unordered maps
			Float_t Erec_E1   = (E1_map.count(EventID_cast) > 0)   ? E1_map[EventID_cast]   : 0;
			Float_t Erec_Eres = (Eres_map.count(EventID_cast) > 0) ? Eres_map[EventID_cast] : 0;

			// Gaussian blur each detector with detector resolution
			Erec_DE_bl = getEres(Erec_DE, 0.8);		
			Erec_E1_bl = getEres(Erec_E1, 1.1);
			Erec_Eres_bl = getEres(Erec_Eres, GAGG_resolution(Erec_Eres));

			// fill banana histgorams
			hDEvE1Ban->Fill(Erec_E1_bl, Erec_DE_bl);
			hDEvEresBan->Fill(Erec_Eres_bl, Erec_DE_bl+Erec_E1_bl);
			// After calculating theta_calc, fill theta specific histograms
			for (int i = 0; i < nThetaBins; ++i) {
			    if (theta_calc >= thetaBinEdges[i] && theta_calc < thetaBinEdges[i+1]) {
			        hDEvE1Ban_thetaBins[i]->Fill(Erec_E1_bl, Erec_DE_bl);		
			        hDEvEresBan_thetaBins[i]->Fill(Erec_Eres_bl, Erec_DE_bl+Erec_E1_bl);
			        break;
			    }
			}
			
			// Calculate total energy by including detector resolution
			if (Erec_DE!=0 && Erec_E1==0 && Erec_Eres==0) { 
				En_tot = Erec_DE; 
				//En_tot += pInSS_interp(En_tot); 	// add back predicted loss in Stainless window
			} else if (Erec_DE!=0 && Erec_E1!=0 && Erec_Eres==0) { 
				En_tot = Erec_DE + Erec_E1; 
				//En_tot += pInSS_interp(En_tot);
			} else if (Erec_DE!=0 && Erec_E1!=0 && Erec_Eres!=0) { 
				En_tot = Erec_DE + Erec_E1 + Erec_Eres; 
				//En_tot += pInSS_interp(En_tot); 
			} else { 
				En_tot = 0; 
				std::cout << "No energy in any  telescope element! Why did this trigger???" << std::endl; 
			}
			
			// Calculating excitation energy
			P_ejc = sqrt(pow((En_tot) + mass_ejectile, 2) - pow(mass_ejectile, 2));
			Exc_en = sqrt(pow(Ek_beam + mass_beam + mass_target - En_tot - mass_ejectile, 2.) - pow(P_beam, 2.) - pow(P_ejc, 2.) + 2*P_beam*P_ejc*cos(theta_calc*M_PI/180)) - mass_recoil; 
	
			// Filling histograms
			hEex->Fill(Exc_en);		hTheta->Fill(theta_calc);		hEprot->Fill(En_tot);

			// Fill accuracy histograms if true values are available
			if (eventMap.count(EventID_cast) > 0) {
				auto true_values = eventMap.at(EventID_cast);
				double theta_true = true_values.first;
				double E_true = true_values.second;
				
				// Fill accuracy histograms (reconstructed - true)
				hThetaAccur->Fill(theta_calc - theta_true);		hThetaAccurvsTheta->Fill(theta_true, theta_calc - theta_true);
				hEprotAccur->Fill(En_tot - E_true);				hEprotAccurvsEprot->Fill(En_tot, En_tot - E_true);
			}

			// check for veto coincidence
			for (int j = 0; j < 4; ++j) {
				if (Erec_DE!=0 && Erec_E1!=0 && E1Veto_list[j].count(EventID_cast)) {
					hE1Veto_dE_posMap->Fill(x_DE, y_DE);
				} else if (Erec_DE!=0 && Erec_E1!=0 && Erec_Eres!=0 && EresVeto_list[j].count(EventID_cast)) {
					hEresVeto_dE_posMap->Fill(x_DE, y_DE);
				}
			}
	
			// check for Dipole and/or HR coincidence
			if (Dip_map.count(EventID_cast)) {
				std::tuple<Float_t,Float_t,Float_t> posTup = Dip_map[EventID_cast];
				x_Dip = std::get<0>(posTup);		y_Dip = std::get<1>(posTup);
				hDipCoinc->Fill(x_Dip, y_Dip);
			}
			
			if (HR_map.count(EventID_cast)) {
				std::tuple<Float_t,Float_t,Float_t> posTup = HR_map[EventID_cast];
				x_HR = std::get<0>(posTup);		y_HR = std::get<1>(posTup);
				hHRcoinc->Fill(x_HR, y_HR);

				// fill banana histgorams
				hHRvEexcCoinc->Fill(x_HR, Exc_en);
				// After calculating theta_calc, fill theta specific histograms
				for (int j = 0; j < nThetaBins; ++j) {
				    if (theta_calc >= thetaBinEdges[j] && theta_calc < thetaBinEdges[j+1]) {
				        hHRvEexcCoinc_thetaBins[j]->Fill(x_HR, Exc_en);
				        break;
				    }
				}
			}
		}
	}
}



// main analysis function
void analysis_nChamb(Int_t excLabel, const char* recType)
{
	// Initialize reaction information
	initializeReactionInfo();
	
	// Use excitation energies from INI file
	if (excLabel >= reactionInfo.recoil_excEns.size()) {
		std::cerr << "Error: excLabel (" << excLabel << ") is out of range for excitation energies array (size: " << reactionInfo.recoil_excEns.size() << ")" << std::endl;
		return;
	}
	Float_t excEn = reactionInfo.recoil_excEns[excLabel];
	
	// _______________________________________________________________________________________________
	// INPUT/OUTPUT FILES, TREES, HISTOGRAMS, VARIABLE DECLARATION ----


	// ========= Input files ==========
  	TFile *ejectile_file = TFile::Open(Form("../%s_sim/Detector_output/Detectors_%s_%s_excEn%02d_ejectile.root", reaction, reaction, recType, excLabel));
	TFile *recoil_file = TFile::Open(Form("../%s_sim/Detector_output/Detectors_%s_%s_excEn%02d_recoil.root", reaction, reaction, recType, excLabel));

	// ========= Define telescope detector structure ==========
	std::vector<TString> detector_names = {"Front", "Back"};		std::vector<TString> detector_prefix = {"h_FP_Eex", "h_BP_Eex"};
	std::vector<TTree*> trees_DE, trees_E1, trees_Eres;
	std::vector<TTree*>  trees_fpE1_veto, trees_bpE1_veto, trees_fpEres_veto, trees_bpEres_veto;
	TTree* dipole_testblock = (TTree*)recoil_file->Get("VirtualDetector/dipole_testblock");
	TTree* HR_testblock = (TTree*)recoil_file->Get("Detector/HR_testblock");
	
	// ========= Define histograms ==========
   	TString hName, hTitle;
	std::vector<TH1D*> hEex, hTheta, hEprot;					// light ejectile properties
	std::vector<TH1D*> hThetaAccur, hEprotAccur;				// accuracy measurement histograms
	std::vector<TH2D*> hThetaAccurvsTheta, hEprotAccurvsEprot;	// 2D accuracy histograms
	std::vector<TH2D*> hDEvE1Ban, hDEvEresBan;					// telescope banana histograms
	std::vector<TH2D*> hE1Veto_dE_posMap, hEresVeto_dE_posMap;	// telescope veto dE position histograms
	std::vector<TH2D*> hDipCoinc, hHRcoinc, hHRvEexcCoinc;		// HR coincidence histograms
	
	for(int i = 0; i < detector_names.size(); ++i) {
		// telescope trees
		trees_DE.push_back(   (TTree*)ejectile_file->Get(Form("Detector/%s_tele_DE", detector_names[i].Data())) );
		trees_E1.push_back(   (TTree*)ejectile_file->Get(Form("Detector/%s_tele_E1", detector_names[i].Data())) );
		trees_Eres.push_back( (TTree*)ejectile_file->Get(Form("Detector/%s_tele_Eres", detector_names[i].Data())) );

		// telescope veto trees
		if(detector_names[i] == "Front") {
			trees_fpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escFP_E1_top") );
			trees_fpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escFP_E1_bot") );
			trees_fpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escFP_E1_inring") );
			trees_fpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escFP_E1_outring") );
			trees_fpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escFP_Eres_top") );
			trees_fpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escFP_Eres_bot") );
			trees_fpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escFP_Eres_inring") );
			trees_fpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escFP_Eres_outring") );
		} else if(detector_names[i] == "Back") {
			trees_bpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escBP_E1_top") );
			trees_bpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escBP_E1_bot") );
			trees_bpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escBP_E1_inring") );
			trees_bpE1_veto.push_back(   (TTree*)ejectile_file->Get("VirtualDetector/escBP_E1_outring") );
			trees_bpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escBP_Eres_top") );
			trees_bpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escBP_Eres_bot") );
			trees_bpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escBP_Eres_inring") );
			trees_bpEres_veto.push_back( (TTree*)ejectile_file->Get("VirtualDetector/escBP_Eres_outring") );
		}

		// light ejectile properties
		hName = labelFromFloat(detector_prefix[i], excEn, "_Eexc_spect", 1);	
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Reconstructed excitation energy; E* (MeV) ; Counts", 1);
		hEex.push_back( new TH1D(hName, hTitle, 100,excEn-5,excEn+5) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_Theta_spect", 1);	
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Reconstructed theta angle; theta (deg) ; Counts", 1);
		hTheta.push_back( new TH1D(hName, hTitle, 180,0,90) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_Eprot_spect", 1);	
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Exact proton energy; Ep (MeV) ; Counts", 1);
		hEprot.push_back( new TH1D(hName, hTitle, 300,0,150) );

		// accuracy histograms
		hName = labelFromFloat(detector_prefix[i], excEn, "_Theta_Accur", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Theta Reconstruction Accuracy; Theta_rec - Theta_true (deg) ; Counts", 1);
		hThetaAccur.push_back( new TH1D(hName, hTitle, 200, -5, 5) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_Theta_Accur_vs_Theta", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Theta Reconstruction Accuracy; Theta (deg); Theta_rec - Theta_true (deg)", 1);
		hThetaAccurvsTheta.push_back( new TH2D(hName, hTitle, 180, 0, 90, 200, -5, 5) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_Eprot_Accur", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Proton Energy Reconstruction Accuracy; Eprot_rec - Eprot_true (MeV) ; Counts", 1);
		hEprotAccur.push_back( new TH1D(hName, hTitle, 200, -5, 5) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_Eprot_Accur_vs_Eprot", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Eprot Reconstruction Accuracy; Eprot (MeV); Eprot_rec - Eprot_true (MeV)", 1);
		hEprotAccurvsEprot.push_back( new TH2D(hName, hTitle, 300, 0, 150, 200, -5, 5) );
		
		// telescope banana histograms
		hName = labelFromFloat(detector_prefix[i], excEn, "_dE_v_E1_Banana", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: dE vs E1 Banana; E1 (MeV) ; dE (MeV)", 1);
		hDEvE1Ban.push_back( new TH2D(hName, hTitle, 500,0,30,500,0,15) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_dE_v_Eres_Banana", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: dE vs Eres Banana; Eres (MeV) ; dE+E1 (MeV)", 1);
		hDEvEresBan.push_back( new TH2D(hName, hTitle, 500,0,150,500,0,50) );

		// telescope veto histograms
		hName = labelFromFloat(detector_prefix[i], excEn, "_E1Veto_dE_posMap", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: dE position map with veto from Escaping E1; x position (mm) ; y position (mm) ; Counts", 1);
		hE1Veto_dE_posMap.push_back( new TH2D(hName, hTitle, 140,-70,70,50,-25,25) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_EresVeto_dE_posMap", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: dE position map with veto from Escaping Eres; x position (mm) ; y position (mm) ; Counts", 1);
		hEresVeto_dE_posMap.push_back( new TH2D(hName, hTitle, 140,-70,70,50,-25,25) );

		// heavy residues coincidences
		hName = labelFromFloat(detector_prefix[i], excEn, "_Dipole_coinc", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: Dipole entrance coincidence; x position (mm) ; y position (mm)", 1);
		hDipCoinc.push_back( new TH2D(hName, hTitle, 400,-200,200,200,-100,100) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_HR_coinc", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: HR plane coincidence; x position (mm) ; y position (mm)", 1);
		hHRcoinc.push_back( new TH2D(hName, hTitle, 400,-200,200,200,-100,100) );
		hName = labelFromFloat(detector_prefix[i], excEn, "_HR_X_v_Eexc_coinc", 1);
		hTitle = labelFromFloat(TString(detector_names[i] + " pocket E* = "), excEn, " MeV: HR plane: X deflection vs Excitation Energy; x position (mm) ; Eexc (MeV)", 1);
		hHRvEexcCoinc.push_back( new TH2D(hName, hTitle, 250,-200,50,300,0,30) );
	}
	std::vector<std::vector<TTree*>> trees_E1_veto = {trees_fpE1_veto, trees_bpE1_veto};
	std::vector<std::vector<TTree*>> trees_Eres_veto = {trees_fpEres_veto, trees_bpEres_veto};

	// histograms as a function of theta
	const Int_t nThetaBins = 9;  	Double_t thetaBinEdges[nThetaBins+1];
	for (int i = 0; i <= nThetaBins; ++i) thetaBinEdges[i] = i*(90/nThetaBins);

	std::vector<TH2D*> hDEvE1Ban_thetaBins;				std::vector<TH2D*> hDEvEresBan_thetaBins;		std::vector<TH2D*> hHRvEexcCoinc_thetaBins;
	for (int i = 0; i < nThetaBins; ++i) {
		TString hname, htitle;

	    hname = labelFromFloat("h_Eex", excEn, Form("_dE_v_E1_Banana_theta_%d_%d", int(thetaBinEdges[i]), int(thetaBinEdges[i+1])), 1);
	    htitle = labelFromFloat("dE v E1 Banana (", excEn, Form(" MeV, %.0f < #theta < %.0f); E1 (MeV); dE (MeV)", thetaBinEdges[i], thetaBinEdges[i+1]), 1);
	    hDEvE1Ban_thetaBins.push_back(new TH2D(hname, htitle, 500, 0, 30, 500, 0, 15));
	    
		hname = labelFromFloat("h_Eex", excEn, Form("_dE_v_ERes_Banana_theta_%d_%d", int(thetaBinEdges[i]), int(thetaBinEdges[i+1])), 1);
	    htitle = labelFromFloat("dE v ERes Banana (E*=", excEn, Form(" MeV, %.0f < #theta < %.0f); Eres (MeV); dE+E1 (MeV)", thetaBinEdges[i], thetaBinEdges[i+1]), 1);
	    hDEvEresBan_thetaBins.push_back(new TH2D(hname, htitle, 500, 0, 150, 500, 0, 50));

		hname = labelFromFloat("hEex", excEn, Form("_HR_X_v_Eexc_coinc_theta_%d_%d", int(thetaBinEdges[i]), int(thetaBinEdges[i+1])), 1);
	    htitle = labelFromFloat("HR plane: X deflection vs Excitation Energy (", excEn, Form(" MeV, %.0f < #theta < %.0f); X deflection (mm); Eexc (MeV)", thetaBinEdges[i], thetaBinEdges[i+1]), 1);
	    hHRvEexcCoinc_thetaBins.push_back(new TH2D(hname, htitle, 250, -200, 50, 300, 0, 30));
	}

	// create writing file. Can also use open_or_create_rootfile function to alternatively just open
	TString fname = Form("../%s_sim/Hist_output/histograms_%s_%s_posFO_targ2.5mm.root", reaction, reaction, recType);
	TFile* file_temp = open_or_create_rootfile(fname.Data());
	
	// ========= Analyse telsecopes ==========
	// Loop over front/back
	for (size_t i = 0; i < detector_names.size(); ++i) 
	{
		analyseTelescope( reaction, recType, excLabel, detector_names[i],
			trees_DE[i],	trees_E1[i],		trees_Eres[i],			trees_E1_veto[i],			trees_Eres_veto[i],		dipole_testblock, 	HR_testblock,
			offsets[i],		rotations[i], 		nThetaBins,				thetaBinEdges,
		  	hEex[i],		hTheta[i],			hEprot[i],				hThetaAccur[i],				hEprotAccur[i],			hThetaAccurvsTheta[i],		hEprotAccurvsEprot[i],
			hDEvE1Ban[i],	hDEvEresBan[i],		hE1Veto_dE_posMap[i],	hEresVeto_dE_posMap[i],		hDEvE1Ban_thetaBins,	hDEvEresBan_thetaBins,
			hDipCoinc[i],	hHRcoinc[i],		hHRvEexcCoinc[i],		hHRvEexcCoinc_thetaBins
	  	);

		hEex[i]->Write();			hTheta[i]->Write();			hEprot[i]->Write();
		hThetaAccur[i]->Write();	hEprotAccur[i]->Write();	hThetaAccurvsTheta[i]->Write();		hEprotAccurvsEprot[i]->Write();
		hDEvE1Ban[i]->Write();		hDEvEresBan[i]->Write();	hE1Veto_dE_posMap[i]->Write();		hEresVeto_dE_posMap[i]->Write();
		hDipCoinc[i]->Write();		hHRcoinc[i]->Write();		hHRvEexcCoinc[i]->Write();
		if(i==detector_names.size()-1) { for (int j = 0; j < nThetaBins; ++j) { 
			hDEvE1Ban_thetaBins[j]->Write();		hDEvEresBan_thetaBins[j]->Write(); 			hHRvEexcCoinc_thetaBins[j]->Write(); 
		} }
	}
	
	file_temp->Close();
}