
// _______________________________________________________________________________________________
// LIBRARIES ----

// ========= C libraries =========
#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <cmath>
#include <math.h>
using namespace std;

// ========= ROOT libraries =========
#include <TROOT.h>
#include <TLatex.h>
#include <TRandom3.h>

// ========= Functions =========
// Adds telescope energy resolution gaussian to energy values
Float_t getEres(Float_t EnergyMeV, Float_t Res_Percent){
	
	Float_t E_res = gRandom->Gaus(EnergyMeV,Res_Percent/100*EnergyMeV);	
	return E_res;
}
// Get mass A and Z from PDGid value
Int_t getZ(Float_t PDGid){

	Int_t Z = floor((PDGid-1000000000)*0.0001);	
	return Z;
}
Int_t getA(Float_t PDGid){

	Int_t A = (PDGid-1000000000-getZ(PDGid)*10000)*0.1;	
	return A;
}

int removeDuplicates(Float_t arr[], Int_t n){ 
 
    // Return, if array is empty or contains a single element 
    if (n == 0 || n == 1) return n; 
    int temp[n];  
    // Start traversing elements 
    int j = 0; 
    // If current element is not equal  to next element then store that  current element 
    for (int i = 0; i < n - 1; i++) 
        if (arr[i] != arr[i + 1]) 
            temp[j++] = arr[i];   
    // Store the last element as whether  it is unique or repeated, it hasn't  stored previously 
    temp[j++] = arr[n - 1];   
    // Modify original array 
    for (int i = 0; i < j; i++) 
        arr[i] = temp[i];   
    return j; 
} 

double getTheta(Float_t x, Float_t y, Float_t z){

	float r=sqrt(pow(-x,2.0)+pow(y,2.0)+pow(z,2.0));
	double theta=acos((double)((z)/r))*180/M_PI;
	return theta;

}


void TreeAnalysis(Float_t mass_beam, Float_t mass_target, Float_t mass_ejectile, Float_t mass_recoil, Float_t EuA_beam){

	// _______________________________________________________________________________________________
	// INPUT/OUTPUT FILES, TREES, HISTOGRAMS, VARIABLE DECLARATION ----

	// ========= Input files ==========
  	// NOTE: Update reaction name in path as needed (e.g., "../206Pb_dd_sim/Detector_output/...")
  	TFile *File_detector = TFile::Open("../Detector_output/Detectors_dd_fis_50000.root"); 
  	//TFile *File_detector = TFile::Open("../Results/238Udp/Detectorsdp.root"); 
	//TFile *SW_deut_file = new TFile("SW_deut_file.root","RECREATE");
  	//TTree *SW_deut_tree = new TTree("SW_deut_tree","Stainless steel window correction for deuterons");
	char name[200];    

  	// ========= Creating trees from file ==========
  	// Telescope
  	TTree *Telescope_DE  = (TTree*)File_detector->Get("Detector/Telescope_DE");
  	TTree *Telescope_E1  = (TTree*)File_detector->Get("Detector/Telescope_E1");
  	TTree *Telescope_E2  = (TTree*)File_detector->Get("Detector/Telescope_E2");
  	TTree *Telescope_E3  = (TTree*)File_detector->Get("Detector/Telescope_E3");
  	TTree *Telescope_E4  = (TTree*)File_detector->Get("Detector/Telescope_E4");
  	TTree *Telescope_E5  = (TTree*)File_detector->Get("Detector/Telescope_E5");
  	TTree *Telescope_E6  = (TTree*)File_detector->Get("Detector/Telescope_E6");
  	TTree *Telescope_E[6] = {Telescope_E1,Telescope_E2,Telescope_E3,Telescope_E4,Telescope_E5,Telescope_E6};
  	//TTree *Telescope_SW  = (TTree*)File_detector->Get("Detector/Telescope_window");

  	// HR detector
  	TTree *HR_Detector = (TTree*)File_detector->Get("Detector/HR_detector");

  	// Fission fragment detectors
  	TTree *Solarcells_top = (TTree*)File_detector->Get("Detector/Solarcells_top");
   	TTree *Silicon_bottom = (TTree*)File_detector->Get("Detector/Silicon_bottom");
   	TTree *Silicon_side = (TTree*)File_detector->Get("Detector/Silicon_side");
   	TTree *FF_Detectors[3] = {Solarcells_top,Silicon_bottom,Silicon_side};


   	// ========= Histograms ==========
   	// Telescope
	TH2D *H_DE_Etot = new TH2D("H_DE_Etot","#DeltaE-E scatter plot; Energy Etot (MeV); Energy #DeltaE (MeV)",512,0.0,50,512,0.0,10.0);
	TH2D *H_DE_Eex_Ekin = new TH2D("H_DE_Eex_Ekin","Excitation energy vs kinetic energy of target-like residues ; E* (MeV) ; Kinetic energy (MeV)",350,0.0,35,550,0,55);

   	TH2D *H_DE_E1 = new TH2D("H_DE_E1","#DeltaE-E1; Energy E1 (MeV); Energy #DeltaE (MeV)",256,0,50,256,0.0,10.0);
   	TH2D *H_DE_E2 = new TH2D("H_DE_E2","#DeltaE-E2; Energy E2 (MeV); Energy #DeltaE (MeV)",256,0,50,256,0.0,10.0);
   	TH2D *H_DE_E3 = new TH2D("H_DE_E3","#DeltaE-E3; Energy E3 (MeV); Energy #DeltaE (MeV)",256,0,50,256,0.0,10.0);
   	TH2D *H_DE_E4 = new TH2D("H_DE_E4","#DeltaE-E4; Energy E4 (MeV); Energy #DeltaE (MeV)",256,0,50,256,0.0,10.0);
   	TH2D *H_DE_E5 = new TH2D("H_DE_E5","#DeltaE-E5; Energy E5 (MeV); Energy #DeltaE (MeV)",256,0,50,256,0.0,10.0);
   	TH2D *H_DE_E6 = new TH2D("H_DE_E6","#DeltaE-E6; Energy E6 (MeV); Energy #DeltaE (MeV)",256,0,50,256,0.0,10.0);
   	TH2D *H_DE_E[6] = {H_DE_E1,H_DE_E2,H_DE_E3,H_DE_E4,H_DE_E5,H_DE_E6};  

   	TH2D *H_E2_E1 = new TH2D("H_E2_E1","E2 vs E1; Energy E2 (MeV); Energy E1 (MeV)",256,0,25,256,0,25);
   	TH2D *H_E3_E2 = new TH2D("H_E3_E2","E3 vs E2; Energy E3 (MeV); Energy E2 (MeV)",256,0,25,256,0,25);
   	TH2D *H_E4_E3 = new TH2D("H_E4_E3","E4 vs E3; Energy E4 (MeV); Energy E3 (MeV)",256,0,25,256,0,25);
   	TH2D *H_E5_E4 = new TH2D("H_E5_E4","E5 vs E4; Energy E5 (MeV); Energy E4 (MeV)",256,0,25,256,0,25);
   	TH2D *H_E6_E5 = new TH2D("H_E6_E5","E6 vs E5; Energy E6 (MeV); Energy E5 (MeV)",256,0,25,256,0,25);
   	TH2D *H_E_E[5] = {H_E2_E1,H_E3_E2,H_E4_E3,H_E5_E4,H_E6_E5};

   	TH2D *H_DE_XY = new TH2D("H_DE_XY","Position of target-like residues; x position (mm); y position (mm)",1600,-10,10,1600,-10,10);
   	TH2D *H_DE_XY_pix = new TH2D("H_DE_XY_pix","Mapping of target-like residues; Vertical strip number; Horizontal strip number",16,1,17,16,1,17);
 	TH1D *H_DE_theta = new TH1D("H_DE_theta","#theta distribution of target-like residues; #theta (deg); Count",512,54,66);
 	TH1D *H_DE_DT = new TH1D("H_DE_DT","#Delta time target-like residues; #Deltat (ns); Count",100,0,5);

 	TH1D *H_Singles_Ekin_VSTRIP_3 = new TH1D("H_Singles_Ekin_VSTRIP_3","Kinetic energy of singles in strip 3; Energy (MeV); Count",512,0,50);
 	TH1D *H_Singles_Eexc_VSTRIP_3 = new TH1D("H_Singles_Eexc_VSTRIP_3","Excitation energy of singles in strip 3; E* (MeV); Count",512,-5,22);

   	// Heavy residues
   	TH1D *H_HR_E = new TH1D("H_HR_E","Kinetic energy of beam-like residues; Energy (MeV); Count",3000,0,3000);
 	TH1D *H_HR_DT = new TH1D("H_HR_DT","#Delta time beam-lke residues; #Deltat (ns); Count",200,350,410);
 	TH1D *H_HR_DT_2 = new TH1D("H_HR_DT_2","#Delta time beam-lke residues; #Deltat (ns); Count",200,350,410);


   	TH2D *H_HR_XY_g = new TH2D("H_HR_XY_g","Mapping of beam-like residues; Vertical strip number; Horizontal strip number",122,0,122,40,-20,20);
    TH2D *H_HR_XY_n = new TH2D("H_HR_XY_n","Mapping of beam-like residues; Vertical strip number; Horizontal strip number",122,0,122,40,-20,20);
   	TH2D *H_HR_XY_2n = new TH2D("H_HR_XY_2n","Mapping of beam-like residues; Vertical strip number; Horizontal strip number",122,0,122,40,-20,20);
   	TH2D *H_HR_XY_3n = new TH2D("H_HR_XY_3n","Mapping of beam-like residues; Vertical strip number; Horizontal strip number",122,0,122,40,-20,20);

   	TH2D *H_HR_Ex_vs_x = new TH2D("H_HR_Ex_vs_x","HR position vs excitation energy; Horizontal position of HR x(mm); CN excitation energy (MeV)",122,0,122,80,0,13);
   
    

  	// Fission fragments
   	TH2D *H_FF_XY_heavy = new TH2D("H_FF_XY_heavy","Fission fragment position;x (mm);y (mm) ",160,-120,120,160,-120,120);
   	TH2D *H_FF_XY_light = new TH2D("H_FF_XY_light","Fission fragment position;x (mm);y (mm) ",160,-120,120,160,-120,120);
 	// Theta angle 1D
	TH1F *HTHETA_top_light = new TH1F("HTHETA_top_light","Theta angle of fission fragments (top det);Theta (deg)",80,0,18);
	TH1F *HTHETA_top_heavy = new TH1F("HTHETA_top_heavy","Theta angle of fission fragments (top det);Theta (deg)",80,0,18);
	TH1F *HTHETA_bot_light = new TH1F("HTHETA_bot_light","Theta angle of fission fragments (top det);Theta (deg)",80,0,18);
	TH1F *HTHETA_bot_heavy = new TH1F("HTHETA_bot_heavy","Theta angle of fission fragments (top det);Theta (deg)",80,0,18);
	TH1F *HTHETA_side_light = new TH1F("HTHETA_side_light","Theta angle of fission fragments (top det);Theta (deg)",80,0,18);
	TH1F *HTHETA_side_heavy = new TH1F("HTHETA_side_heavy","Theta angle of fission fragments (top det);Theta (deg)",80,0,18);
	// Theta angle 1D ve energy 2D
	TH2F *HThetaE_top_light = new TH2F("HThetaE_top_light","Energy per angle;Theta (deg);E (MeV)",40,0,20,256,0,3000);
	TH2F *HThetaE_top_heavy = new TH2F("HThetaE_top_heavy","Energy per angle;Theta (deg);E (MeV)",40,0,20,256,0,3000);
	TH2F *HThetaE_bot_light = new TH2F("HThetaE_bot_light","Energy per angle;Theta (deg);E (MeV)",40,0,20,256,0,3000);
	TH2F *HThetaE_bot_heavy = new TH2F("HThetaE_bot_heavy","Energy per angle;Theta (deg);E (MeV)",40,0,20,256,0,3000);
	TH2F *HThetaE_side_light = new TH2F("HThetaE_side_light","Energy per angle;Theta (deg);Fission fragment energy (MeV)",40,0,20,256,0,3000);
	TH2F *HThetaE_side_heavy = new TH2F("HThetaE_side_heavy","Energy per angle;Theta (deg);Fission fragment energy (MeV)",40,0,20,256,0,3000);
	// Delta T 1D
	TH1D *H_FF_DT_top = new TH1D("H_FF_DT_top","#Delta time fission fragments;dt (ns)",100,0,14);
	TH1D *H_FF_DT_bot = new TH1D("H_FF_DT_bot","#Delta time fission fragments;dt (ns)",100,0,14);
	TH1D *H_FF_DT_side = new TH1D("H_FF_DT_side","#Delta time fission fragments;dt (ns)",100,0,14);
	// Theta angle vs position 2D
	//TH2F *HyE_top_light = new TH2F("HyE_top_light","Energy per y position top detector;y (mm);E (MeV)",16,24,64,2500,0,2500);
	//TH2F *HyE_top_heavy = new TH2F("HyE_top_heavy","Energy per y position top detector;y (mm);E (MeV)",16,24,64,2500,0,2500);
	//TH2F *HyE_bot_light = new TH2F("HyE_bot_light","Energy per y position bottom detector;y (mm);E (MeV)",16,-24,-64,2500,0,2500);
	//TH2F *HyE_bot_heavy = new TH2F("HyE_bot_heavy","Energy per y position bottom detector;y (mm);E (MeV)",8,-24,-64,2500,0,2500);
	//TH2F *HxE_side_light = new TH2F("HxE_side_light","Energy per x position side detector;x (mm);E (MeV)",80,0,-100,2500,0,2500);
	//TH2F *HxE_side_heavy = new TH2F("HxE_side_heavy","Energy per x position side detector;x (mm);E (MeV)",50,0,-100,2500,0,2500);

	TH1D *H_FF_E_top = new TH1D("H_FF_E_top","Energy of fission fragments; Energy (MeV)",100,0,3000);
	TH1D *H_FF_E_bot = new TH1D("H_FF_E_bot","Energy of fission fragments; Energy (MeV)",100,0,3000);
	TH1D *H_FF_E_side = new TH1D("H_FF_E_side","Energy of fission fragments; Energy (MeV)",100,0,3000);

	TH2F *HyE_top_light = new TH2F("HyE_top_light","Energy per y position top detector;H strip n; E (MeV)",16,0.5,16.5,256,0,3000);
	TH2F *HyE_top_heavy = new TH2F("HyE_top_heavy","Energy per y position top detector;H strip n; E (MeV)",16,0.5,16.5,256,0,3000);
	TH2F *HyE_bot_light = new TH2F("HyE_bot_light","Energy per y position bottom detector;V strip n; E (MeV)",16,0.5,16.5,256,0,3000);
	TH2F *HyE_bot_heavy = new TH2F("HyE_bot_heavy","Energy per y position bottom detector;H strip n; E (MeV)",16,0.5,16.5,256,0,3000);
	TH2F *HxE_side_light = new TH2F("HxE_side_light","Energy per x position side detector;Side fission detector vertical strip (n); Fission fragment energy (MeV)",122,0.5,122.5,256,0,3000);
	TH2F *HxE_side_heavy = new TH2F("HxE_side_heavy","Energy per x position side detector;Side fission detector vertical strip (n); Fission fragment energy (MeV)",122,0.5,122.5,256,0,3000);
 
   	TH2D *H_FF_XY_heavy_top = new TH2D("H_FF_XY_heavy_top","Fission fragment position top detector;V strip n;V strip n ",16,0,16,16,0,16);
   	TH2D *H_FF_XY_light_top = new TH2D("H_FF_XY_light_top","Fission fragment position top detector;V strip n;V strip n",16,0,16,16,0,16);		
   	TH2D *H_FF_XY_heavy_bot = new TH2D("H_FF_XY_heavy_bot","Fission fragment position bot detector;V strip n;V strip n",16,0,16,16,0,16);
   	TH2D *H_FF_XY_light_bot = new TH2D("H_FF_XY_light_bot","Fission fragment position bot detector;V strip n;V strip n",16,0,16,16,0,16);
   	TH2D *H_FF_XY_heavy_side = new TH2D("H_FF_XY_heavy_side","Fission fragment position side detector;V strip n;V strip n",122,0,122.5,40,0.5,40.5);
   	TH2D *H_FF_XY_light_side = new TH2D("H_FF_XY_light_side","Fission fragment position side detector;V strip n;V strip n",122,0,122.5,40,0.5,40.5);

    TH2F *h2_HRmap_Ex_strip[16];
    for(int i=0;i<16;i++){
        sprintf(name,"h2_HRmap_Ex_strip%d",i+1);
        h2_HRmap_Ex_strip[i] = new TH2F(name,name,120,0,120,75,-3,15);        
    }


  	// ========= Variable declaration ==========
 	// Beam
   	Float_t Ek_beam = mass_beam/931.5*EuA_beam;
	Float_t P_beam = sqrt(pow(Ek_beam+mass_beam,2)-pow(mass_beam,2));  	
   	// Telescope
   	Int_t Tel_VSTRIP, Tel_HSTRIP;
   	Float_t x_DE_pix, y_DE_pix, z_DE_pix=99.025,theta_DE_pix;
	Float_t Etot_E=0, Energy_Ex=0, Etot_tel_res, P_tel_res, Exc_energy, SW_energy=0;   		
	// Heavy residue
	Int_t A_HR;
	// Fission



   	gStyle->SetPalette(109);	// Set palette style DeepSea
    UInt_t SigStep=0;

   	// _______________________________________________________________________________________________
	// BRANCHING TREES ----

   	// ========= Variables for tree branches ==========

   	// Telescope - DeltaE detector
   	Long64_t nEntries_DE;
   	Float_t x_DE, y_DE, z_DE, Px_DE, Py_DE, Pz_DE, t_DE, PDGid_DE, EventID_DE, Edep_DE;
   	// Telescope - E detector stack
   	Float_t nEntries_E[6], t_E[6], PDGid_E[6], EventID_E[6],Edep_E[6];
   	// Telescope window
   	Float_t nEntries_SW, Edep_SW, EventID_SW;


   	// Heavy residue 
   	Float_t nEntries_HR, x_HR, y_HR, z_HR, Px_HR, Py_HR, Pz_HR, t_HR, PDGid_HR, EventID_HR, Edep_HR;
   	// Fission fragments
   	Float_t nEntries_FF[3], x_FF[3], y_FF[3], z_FF[3], t_FF[3], PDGid_FF[3], EventID_FF[3], Edep_FF[3];
   	Int_t FFtop_VSTRIP, FFtop_HSTRIP;
   	Int_t FFbot_VSTRIP, FFbot_HSTRIP;
   	Int_t FFside_VSTRIP, FFside_HSTRIP;

   	// ========= Tree branching ==========

  	// Telescope - DeltaE detector branching
  	nEntries_DE = Telescope_DE->GetEntries(); 
  	Telescope_DE->SetBranchAddress("x",&x_DE);		Telescope_DE->SetBranchAddress("y",&y_DE);		Telescope_DE->SetBranchAddress("z",&z_DE);
  	Telescope_DE->SetBranchAddress("Px",&Px_DE);	Telescope_DE->SetBranchAddress("Py",&Py_DE);	Telescope_DE->SetBranchAddress("Pz",&Pz_DE);
  	Telescope_DE->SetBranchAddress("t",&t_DE);		Telescope_DE->SetBranchAddress("Edep",&Edep_DE);
  	Telescope_DE->SetBranchAddress("PDGid",&PDGid_DE);	Telescope_DE->SetBranchAddress("EventID",&EventID_DE);
 	
 	// Telescope - E detector stack branching
   	for(Int_t i=0;i<6;i++){

   		nEntries_E[i] = Telescope_E[i]->GetEntries();
   		Telescope_E[i]->SetBranchAddress("t",&t_E[i]);
   		Telescope_E[i]->SetBranchAddress("Edep",&Edep_E[i]);
   		Telescope_E[i]->SetBranchAddress("PDGid",&PDGid_E[i]);
   		Telescope_E[i]->SetBranchAddress("EventID",&EventID_E[i]);
   		
   	}
/*
  	// Telescope window
  	nEntries_SW = Telescope_SW->GetEntries(); 
  	Telescope_SW->SetBranchAddress("Edep",&Edep_SW);
  	Telescope_SW->SetBranchAddress("EventID",&EventID_SW);
  	SW_deut_tree->Branch("Etot_tel_res",&Etot_tel_res,"Etot_tel_res");
	SW_deut_tree->Branch("theta_DE_pix",&theta_DE_pix,"theta_DE_pix");
   	SW_deut_tree->Branch("SW_energy",&SW_energy,"SW_energy");		
*/
   	// Heavy residue detector branching
  	nEntries_HR = HR_Detector->GetEntries(); 
  	HR_Detector->SetBranchAddress("x",&x_HR);	HR_Detector->SetBranchAddress("y",&y_HR);	HR_Detector->SetBranchAddress("z",&z_HR);
  	HR_Detector->SetBranchAddress("Px",&Px_HR);	HR_Detector->SetBranchAddress("Py",&Py_HR);	HR_Detector->SetBranchAddress("Pz",&Pz_HR);
  	HR_Detector->SetBranchAddress("t",&t_HR);	HR_Detector->SetBranchAddress("Edep",&Edep_HR);
  	HR_Detector->SetBranchAddress("PDGid",&PDGid_HR);	HR_Detector->SetBranchAddress("EventID",&EventID_HR);

 	// FF detectors branching
   	for(Int_t i=0;i<3;i++){

   		nEntries_FF[i] = FF_Detectors[i]->GetEntries();
   		FF_Detectors[i]->SetBranchAddress("x",&x_FF[i]);	FF_Detectors[i]->SetBranchAddress("y",&y_FF[i]);	FF_Detectors[i]->SetBranchAddress("z",&z_FF[i]);
   		FF_Detectors[i]->SetBranchAddress("t",&t_FF[i]);
   		FF_Detectors[i]->SetBranchAddress("Edep",&Edep_FF[i]);
   		FF_Detectors[i]->SetBranchAddress("PDGid",&PDGid_FF[i]);
   		FF_Detectors[i]->SetBranchAddress("EventID",&EventID_FF[i]);
   		
   	}

   	// _______________________________________________________________________________________________
	// TRANSFORMING FISSION DOUBLE EVENTS INTO SINGLE EVENTS ----
  
   	Float_t FFevents_sum = nEntries_FF[0]+nEntries_FF[1]+nEntries_FF[2];
   	Float_t EventsFF_ID_raw[(Int_t)FFevents_sum];   	
   	Int_t maxID=0;
   	Float_t EventsFF_ID_single[(Int_t)FFevents_sum];
   	Float_t EventsFF_ID_double[(Int_t)FFevents_sum];
   	Int_t RandHeavyLight;
   	Int_t i1=0, i2=0, i3=0, ndouble = 0, nsingle = 0;	

	// Filling an array of all eventIDs
  	for(Int_t i=0;i<3;i++){ 		
  		for(Int_t nevent_FF=0;nevent_FF<nEntries_FF[i];nevent_FF++){
  			FF_Detectors[i]->GetEntry(nevent_FF);			
  			EventsFF_ID_raw[maxID+nevent_FF]=EventID_FF[i];
  		}
  		maxID=maxID+nEntries_FF[i];
  	}

  	// Looping through the array to find events that are once or twice
	for(Int_t Raw_ID=0;Raw_ID<FFevents_sum;Raw_ID++){
		Int_t filled=0;
		for(Int_t Raw_ID2=0;Raw_ID2<FFevents_sum;Raw_ID2++){
			if((EventsFF_ID_raw[Raw_ID]==EventsFF_ID_raw[Raw_ID2]) && (Raw_ID!=Raw_ID2)){
				EventsFF_ID_double[ndouble]=EventsFF_ID_raw[Raw_ID];
				ndouble++;
				filled=1;
				continue;
			}				
		}
		if(filled==0){
			EventsFF_ID_single[nsingle]=EventsFF_ID_raw[Raw_ID];
			//cout<<"single event"<<EventsFF_ID_raw[Raw_ID]<<endl;
			nsingle++;
		}
	}

	// Sorting the arrays of single and double fission event IDs
	sort(EventsFF_ID_double,EventsFF_ID_double + ndouble);
	sort(EventsFF_ID_single,EventsFF_ID_single + nsingle);
	// Removing duplicates from double fission event ID array
	int newsize = removeDuplicates(EventsFF_ID_double,ndouble);

   	// _______________________________________________________________________________________________
	// LOOPING FOR TELESCOPE ----

  	for(Long64_t nevent_DE=0;nevent_DE<nEntries_DE;nevent_DE++){

        if(nevent_DE%(nEntries_DE/100)==0 && (nevent_DE-1)%(nEntries_DE/100)!=0) {
            cout<<"~~~~~ Progress : "<< SigStep++<<"%\r";
            fflush(stdout);
        }

  		Telescope_DE->GetEntry(nevent_DE);
  		//cout<<"event id "<<EventID_DE;

  		// 1 - Telescope strip calculation from ejectile coordinates
  		Tel_VSTRIP = 17-ceil((x_DE+10)/1.25);		Tel_HSTRIP = ceil((y_DE+10)/1.25);

  		// 2 - Telescope angle calculation pixel coordinate (taken at center of each pixel)
		x_DE_pix = -10+1.25*(Tel_VSTRIP-1)+1.25/2;	y_DE_pix = -10+1.25*(Tel_HSTRIP-1)+1.25/2;
		Float_t alpha1 = (Float_t)(atan(x_DE_pix/z_DE_pix))+60*M_PI/180;
		Float_t rho1 = (Float_t)(sqrt(pow(x_DE_pix,2)+pow(z_DE_pix,2)));	
		Float_t z = (Float_t)(cos(alpha1))*rho1;
		Float_t rho = (Float_t)(sqrt(pow(x_DE_pix,2)+pow(y_DE_pix,2)+pow(z_DE_pix,2)));
		theta_DE_pix = (Float_t)acos(z/rho);

		// 3 - Looping through each E detector to get energy of the E stack

		// Looping through each entry of each E detector

		// Detector E1
		for(Int_t nevent_E1=0;nevent_E1<nEntries_E[0];nevent_E1++)
		{
			Telescope_E[0]->GetEntry(nevent_E1);

			if(EventID_E[0]==EventID_DE)
			{
				Etot_E+=getEres(Edep_E[0],1.3);
				if (Tel_VSTRIP==1) H_DE_E[0]->Fill(getEres(Edep_E[0],1.3),getEres(Edep_DE,1.3));

				// Detector E2
				for(Int_t nevent_E2=0;nevent_E2<nEntries_E[1];nevent_E2++)
				{
					Telescope_E[1]->GetEntry(nevent_E2);

					if(EventID_E[1]==EventID_DE)
					{
						Etot_E+=getEres(Edep_E[1],1.3);
						if (Tel_VSTRIP==1) H_DE_E[1]->Fill(getEres(Edep_E[1],1.3),getEres(Edep_DE,1.3));
						H_E_E[0]->Fill(getEres(Edep_E[1],1.3),getEres(Edep_E[0],1.3));

						// Detector E3
						for(Int_t nevent_E3=0;nevent_E3<nEntries_E[2];nevent_E3++)
						{
							Telescope_E[2]->GetEntry(nevent_E3);

							if(EventID_E[2]==EventID_DE)
							{
								Etot_E+=getEres(Edep_E[2],1.3);
								if (Tel_VSTRIP==1) H_DE_E[2]->Fill(getEres(Edep_E[2],1.3),getEres(Edep_DE,1.3));
								H_E_E[1]->Fill(getEres(Edep_E[2],1.3),getEres(Edep_E[1],1.3));

								// Detector E4
								for(Int_t nevent_E4=0;nevent_E4<nEntries_E[3];nevent_E4++)
								{
									Telescope_E[3]->GetEntry(nevent_E4);

									if(EventID_E[3]==EventID_DE)
									{
										Etot_E+=getEres(Edep_E[3],1.3);
										if (Tel_VSTRIP==1) H_DE_E[3]->Fill(getEres(Edep_E[3],1.3),getEres(Edep_DE,1.3));
										H_E_E[2]->Fill(getEres(Edep_E[3],1.3),getEres(Edep_E[2],1.3));

										for(Int_t nevent_E5=0;nevent_E5<nEntries_E[4];nevent_E5++)
										{
											Telescope_E[4]->GetEntry(nevent_E5);

											if(EventID_E[4]==EventID_DE)
											{
												Etot_E+=getEres(Edep_E[4],1.3);
												if (Tel_VSTRIP==1) H_DE_E[4]->Fill(getEres(Edep_E[4],1.3),getEres(Edep_DE,1.3));
												H_E_E[3]->Fill(getEres(Edep_E[4],1.3),getEres(Edep_E[3],1.3));


												for(Int_t nevent_E6=0;nevent_E6<nEntries_E[5];nevent_E6++)
												{
													Telescope_E[5]->GetEntry(nevent_E6);

													if(EventID_E[5]==EventID_DE)
													{
														Etot_E+=getEres(Edep_E[5],1.3);
														if (Tel_VSTRIP==1) H_DE_E[5]->Fill(getEres(Edep_E[5],1.3),getEres(Edep_DE,1.3));
														H_E_E[4]->Fill(getEres(Edep_E[5],1.3),getEres(Edep_E[4],1.3));

													}
												}
											}
										}
									}
								}
							}	
						}
					}
				}
			}
		}

		// 4 - Calculating excitation energy
		Etot_tel_res = Etot_E+getEres(Edep_DE,1.6);
		P_tel_res = sqrt(pow((Etot_tel_res-Etot_tel_res*0.07)+mass_ejectile,2)-pow(mass_ejectile,2));
		Exc_energy = sqrt(pow(Ek_beam+mass_beam+mass_target-Etot_tel_res-mass_ejectile,2.)-pow(P_beam,2.)-pow(P_tel_res,2.)+2*P_beam*P_tel_res*TMath::Cos(theta_DE_pix))-mass_recoil; 
		/*
		// 5 - Filling output tree with telescope window energy
		for(Int_t nevent_SW=0;nevent_SW<nEntries_SW;nevent_SW++){
			Telescope_SW->GetEntry(nevent_SW);
			if(EventID_SW==EventID_DE){
				SW_energy=Edep_SW;
				SW_deut_tree->Fill();
			}	
		}
		*/
		//cout<<EventID_DE<<endl;
  		// 6 - Filling DeltaE histograms
  		if (!(Tel_VSTRIP==1)&&!(Tel_VSTRIP==16)&&!(Tel_HSTRIP==1)&&!(Tel_HSTRIP==16)) H_DE_Etot->Fill(Etot_E,getEres(Edep_DE,1.6));
  		H_DE_XY_pix->Fill(Tel_VSTRIP,Tel_HSTRIP);

  		
  		H_DE_theta->Fill(theta_DE_pix*180/M_PI);

  		H_DE_DT->Fill(t_DE);
  		H_DE_Eex_Ekin->Fill(Exc_energy,Etot_tel_res);
  		//H_Singles_Eexc_VSTRIP_3->Fill(Exc_energy);

  		if (Tel_VSTRIP==3){

  			H_Singles_Ekin_VSTRIP_3->Fill(Etot_E+Edep_DE);
  			

  		}

  		

  		// _______________________________________________________________________________________________
		// LOOPING FOR HEAVY RESIDUE ----
  		/*
  		// 7 - Looping through each entry of the HR detector
  		for(Int_t nevent_HR=0;nevent_HR<nEntries_HR;nevent_HR++){

  			HR_Detector->GetEntry(nevent_HR);

  			if(EventID_HR==EventID_DE){
  				A_HR = getA(PDGid_HR);
  				//cout<<PDGid_HR<<endl;
  				//printf (" pdgid: %11f \n", PDGid_HR);
  				H_HR_XY_g->Fill(61-x_HR,y_HR);
  				
  				h2_HRmap_Ex_strip[Tel_VSTRIP]->Fill(61-x_HR,Exc_energy);

  				}
  					



  				//H_HR_DT->Fill(t_HR-t_DE);
  				H_HR_E->Fill(Edep_HR);
  				
  				
  		
  		}*/
		
  	  	// _______________________________________________________________________________________________
		// LOOPING FOR FISSION FRAGMENTS ----
  		
  		// 8 - Filling events for top detector
  		for(Int_t nevent_FF0=0;nevent_FF0<nEntries_FF[0];nevent_FF0++){
  			//cout<<"De id "<<EventID_DE<<endl;
  			FF_Detectors[0]->GetEntry(nevent_FF0);
  			FFtop_VSTRIP = 17-ceil((x_FF[0]+40)/5);		FFtop_HSTRIP = ceil((y_FF[0]+20)/2.5);
  			//cout<<nevent_FF0<<endl;

  			for(Int_t j=0;j<nsingle;j++){  				
  				if((EventID_FF[0]==EventID_DE) && (EventID_FF[0]==EventsFF_ID_single[j])){
  					
  					if (PDGid_FF[0]<1000450000){
  						H_FF_E_top->Fill(Edep_FF[0]);
  						//cout<<PDGid_FF[0]<<endl;
  						H_FF_XY_light->Fill(-x_FF[0],y_FF[0]+48);
  						HTHETA_top_light->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1));
  						HThetaE_top_light->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1),getEres(Edep_FF[0],4.0));
  						H_FF_DT_top->Fill(t_FF[0]);
  						HyE_top_light->Fill(FFtop_HSTRIP,getEres(Edep_FF[0],4.0));
  						H_FF_XY_light_top->Fill(FFtop_VSTRIP-0.5,FFtop_HSTRIP-0.5);

  					}
  					else {
  						H_FF_E_top->Fill(Edep_FF[0]);
  						H_FF_XY_heavy->Fill(-x_FF[0],y_FF[0]+48);
  						HTHETA_top_heavy->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1));
  						HThetaE_top_heavy->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1),getEres(Edep_FF[0],4.0));
  						H_FF_DT_top->Fill(t_FF[0]);  
  						HyE_top_heavy->Fill(FFtop_HSTRIP,getEres(Edep_FF[0],4.0));	
  						H_FF_XY_heavy_top->Fill(FFtop_VSTRIP-0.5,FFtop_HSTRIP-0.5);				
  					}
  					continue;
  				}				
			}
	
  			for(Int_t k=0;k<newsize;k++){
  				if((EventID_FF[0]==EventID_DE) && (EventID_FF[0]==EventsFF_ID_double[k]) && ((Int_t)EventID_FF[0]%2==0)){
  					
  					if (PDGid_FF[0]<1000450000){
  						H_FF_E_top->Fill(Edep_FF[0]);
  						H_FF_XY_light->Fill(-x_FF[0],y_FF[0]+48);
  						HTHETA_top_light->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1));
  						HThetaE_top_light->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1),getEres(Edep_FF[0],4.0));
  						H_FF_DT_top->Fill(t_FF[0]);
  						HyE_top_light->Fill(FFtop_HSTRIP,getEres(Edep_FF[0],4.0));
  						H_FF_XY_light_top->Fill(FFtop_VSTRIP-0.5,FFtop_HSTRIP-0.5);
  					}
  					else{
  						H_FF_E_top->Fill(Edep_FF[0]);
  						H_FF_XY_heavy->Fill(-x_FF[0],y_FF[0]+48);
  						HTHETA_top_heavy->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1));
    					HThetaE_top_heavy->Fill(getTheta(-x_FF[0],y_FF[0]+48,z_FF[0]+389.1),getEres(Edep_FF[0],4.0));
  						H_FF_DT_top->Fill(t_FF[0]);
  						HyE_top_heavy->Fill(FFtop_HSTRIP,getEres(Edep_FF[0],4.0));
  						H_FF_XY_heavy_top->Fill(FFtop_VSTRIP-0.5,FFtop_HSTRIP-0.5);
  					}	
  					continue;
  				}				
			}
  		}
		
		
  		// 9 - Filling events for bottom detector
  		for(Int_t nevent_FF1=0;nevent_FF1<nEntries_FF[1];nevent_FF1++){
  			//cout<<"De id "<<EventID_DE<<endl;
  			FF_Detectors[1]->GetEntry(nevent_FF1);
  			FFbot_VSTRIP = 17-ceil((x_FF[1]+40)/5);		FFbot_HSTRIP = 17-ceil((y_FF[1]+20)/2.5);

  			for(Int_t j=0;j<nsingle;j++){
  				if((EventID_FF[1]==EventID_DE) && (EventID_FF[1]==EventsFF_ID_single[j])){
  					
  					if (PDGid_FF[1]<1000450000){
  						H_FF_E_bot->Fill(Edep_FF[1]);
  						H_FF_XY_light->Fill(-x_FF[1],y_FF[1]-47);
  						//cout<<getTheta(-x_FF[1],y_FF[1]-44,z_FF[0]+209.1)<<endl;
  						HTHETA_bot_light->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1));
  						HThetaE_bot_light->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1),getEres(Edep_FF[1],4.0));
  						H_FF_DT_bot->Fill(t_FF[1]);
  						HyE_bot_light->Fill(FFbot_HSTRIP,getEres(Edep_FF[1],4.0));
  						H_FF_XY_light_bot->Fill(FFbot_VSTRIP-0.5,16-FFbot_HSTRIP+1-0.5);
  					}
  					else {
  						H_FF_E_bot->Fill(Edep_FF[1]);
  						H_FF_XY_heavy->Fill(-x_FF[1],y_FF[1]-47);
  						HTHETA_bot_heavy->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1));
  						HThetaE_bot_heavy->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1),getEres(Edep_FF[1],4.0));
  						H_FF_DT_bot->Fill(t_FF[1]);
  						HyE_bot_heavy->Fill(FFbot_HSTRIP,getEres(Edep_FF[1],4.0));
  						H_FF_XY_heavy_bot->Fill(FFbot_VSTRIP-0.5,16-FFbot_HSTRIP+1-0.5);
  					}
  					continue;
  				}				
			}
  					
  			for(Int_t k=0;k<newsize;k++){
  				if((EventID_FF[1]==EventID_DE) && (EventID_FF[1]==EventsFF_ID_double[k]) && ((Int_t)EventID_FF[1]%2!=0)){
  					
  					if (PDGid_FF[1]<1000450000){
  						H_FF_E_bot->Fill(Edep_FF[1]);
  						H_FF_XY_light->Fill(-x_FF[1],y_FF[1]-47);
  						HTHETA_bot_light->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1));
  						HThetaE_bot_light->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1),getEres(Edep_FF[1],4.0));  						
  						H_FF_DT_bot->Fill(t_FF[1]);
  						HyE_bot_light->Fill(FFbot_HSTRIP,getEres(Edep_FF[1],4.0));
  						H_FF_XY_light_bot->Fill(FFbot_VSTRIP-0.5,16-FFbot_HSTRIP+1-0.5);
  					}
  					else {
  						H_FF_E_bot->Fill(Edep_FF[1]);
  						H_FF_XY_heavy->Fill(-x_FF[1],y_FF[1]-46);
  						HTHETA_bot_heavy->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1));
  						HThetaE_bot_heavy->Fill(getTheta(-x_FF[1],y_FF[1]-47,z_FF[1]+209.1),getEres(Edep_FF[1],4.0));  						
  						H_FF_DT_bot->Fill(t_FF[1]);
  						HyE_bot_heavy->Fill(FFbot_HSTRIP,getEres(Edep_FF[1],4.0));
  						H_FF_XY_heavy_bot->Fill(FFbot_VSTRIP-0.5,16-FFbot_HSTRIP+1-0.5);
  					}
  					continue;
  				}				
			}
  		}
		
		
  		// 10 - Filling events for side detector
  		for(Int_t nevent_FF2=0;nevent_FF2<nEntries_FF[2];nevent_FF2++){
  			//cout<<"De id "<<EventID_DE<<endl;
  			FF_Detectors[2]->GetEntry(nevent_FF2);
  			FFside_VSTRIP = 122-ceil((x_FF[2]+61)/1);		FFside_HSTRIP = 20-ceil((y_FF[2])/1);

  			for(Int_t j=0;j<nsingle;j++){
  				if((EventID_FF[2]==EventID_DE) && (EventID_FF[2]==EventsFF_ID_single[j])){
  					
  					 if (PDGid_FF[2]<1000450000){
  					 	H_FF_E_side->Fill(Edep_FF[2]);
  						H_FF_XY_light->Fill(-x_FF[2]-79,y_FF[2]);
  						HTHETA_side_light->Fill(getTheta(-x_FF[2]-79,y_FF[2],z_FF[2]+329.1));
  						HThetaE_side_light->Fill(getTheta(-x_FF[2]-79,y_FF[2],z_FF[2]+329.1),getEres(Edep_FF[2],4.0));
  						H_FF_DT_side->Fill(t_FF[2]);
  						HxE_side_light->Fill(120-FFside_VSTRIP+1,getEres(Edep_FF[2],4.0));
  						H_FF_XY_light_side->Fill(FFside_VSTRIP,FFside_HSTRIP);
  					}
  					else {
  						H_FF_E_side->Fill(Edep_FF[2]);
  						H_FF_XY_heavy->Fill(-x_FF[2]-79,y_FF[2]);
  						HTHETA_side_heavy->Fill(getTheta(-x_FF[2]-79,y_FF[2],z_FF[2]+329.1));
  						HThetaE_side_heavy->Fill(getTheta(-x_FF[2]-79,y_FF[2],z_FF[2]+329.1),getEres(Edep_FF[2],4.0));
  						H_FF_DT_side->Fill(t_FF[2]);
  						HxE_side_heavy->Fill(120-FFside_VSTRIP+1,getEres(Edep_FF[2],4.0));
  						H_FF_XY_heavy_side->Fill(FFside_VSTRIP,FFside_HSTRIP);
  					}
  					continue;
  				}				
			}
  		}
  		

		Etot_E=0;

  	} 

  	// _______________________________________________________________________________________________
	// FISSION DETECION EFFICIENCY ---- 

  	Float_t detected=H_FF_XY_light->Integral()+H_FF_XY_heavy->Integral();
	Float_t detected_top=(double)(H_FF_XY_light_top->Integral()+H_FF_XY_heavy_top->Integral());
	Float_t detected_bot=(double)(H_FF_XY_light_bot->Integral()+H_FF_XY_heavy_bot->Integral());
	Float_t detected_topbot=detected_top+detected_bot;
	Float_t detected_side=(double)(H_FF_XY_light_side->Integral()+H_FF_XY_heavy_side->Integral());
	Float_t not_detected=(nEntries_DE-detected);

/*
	Float_t efficiencies[4] = {detected_top,detected_bot,detected_side,not_detected};
	Int_t piecolors[] = {2,3,4,5};
    TPie *Fission_efficiency = new TPie("Fission_efficiency","Fission event detection efficiency",4,efficiencies,piecolors);
    Fission_efficiency->SetEntryLabel(0,"Top detector");
    Fission_efficiency->SetEntryLabel(1,"Bot detector");
    Fission_efficiency->SetEntryLabel(2,"Side detector");
    Fission_efficiency->SetEntryLabel(3,"Not detected");
*/
	Float_t efficiencies[3] = {detected_topbot,detected_side,not_detected};
	Int_t piecolors[] = {30,40,18};
    TPie *Fission_efficiency = new TPie("Fission_efficiency","Fission detection efficiency",3,efficiencies,piecolors);
    Fission_efficiency->SetEntryLabel(0,"Top+Bot detectors");
    Fission_efficiency->SetEntryLabel(1,"Side detector");
    Fission_efficiency->SetEntryLabel(2,"Not detected");


  	// _______________________________________________________________________________________________
	// CANVAS ----
/*
  	TCanvas *c_de_etot = new TCanvas("c_de_etot","Scatter plots DeltaE-Etot",1000,800);
  	c_de_etot->Divide(2,1);
    c_de_etot->cd(1);
  	H_DE_Etot->Draw("colz");
    c_de_etot->cd(2);
  	H_DE_Eex_Ekin->Draw("colz");


    TCanvas *c_de_ex = new TCanvas("c_de_ex","Scatter plots DeltaE-Ei",1000,800);
    c_de_ex->Divide(3,2);

    for(Int_t j=0;j<6;j++){

    	c_de_ex->cd(j+1);
    	H_DE_E[j]->Draw("colz");
    }


    TCanvas *c_ex_ey = new TCanvas("c_ex_ey","Scatter plots Ei vs Ej",1000,800);
    c_ex_ey->Divide(3,2);

    for(Int_t j=0;j<5;j++){

    	c_ex_ey->cd(j+1);
    	H_E_E[j]->Draw("colz");
    }


    TCanvas *c_de_posangle = new TCanvas("c_de_posangle","Positions and angles of target-like residues",1000,800);
    c_de_posangle->Divide(2,2);
    c_de_posangle->cd(1);   
    H_DE_XY_pix->Draw("colz");
    c_de_posangle->cd(3); 
    H_DE_XY_pix->Draw("LEGO2Z");
    c_de_posangle->cd(2);
    H_DE_theta->Draw();
    c_de_posangle->cd(4);
    H_DE_DT->Draw();

    TCanvas *c_hr_posangle = new TCanvas("c_hr_posangle","Positions and angles of beam-like residues",1200,800);
    c_hr_posangle->Divide(2,2);
    c_hr_posangle->cd(1);
    H_HR_XY_g->Draw("colz");
    c_hr_posangle->cd(3);  
    H_HR_XY_g->Draw("LEGO2Z");  
    c_hr_posangle->cd(2);  
    H_HR_E->Draw();  
    c_hr_posangle->cd(4);     
    H_HR_DT_2->SetLineColor(kBlack);
    H_HR_DT_2->Draw();  
    H_HR_DT->SetLineColor(kGreen);
    H_HR_DT->Draw("same");
*/

    TCanvas *c_FF_pos = new TCanvas("c_FF_pos","Positions and angles of fission fragments",1200,800);
    c_FF_pos->Divide(2,2);
    c_FF_pos->cd(1);
    H_FF_XY_light->Draw("colz");
    H_FF_XY_heavy->Draw("col same");
    c_FF_pos->cd(3);
    H_FF_XY_light->Draw("LEGO2Z");
    H_FF_XY_heavy->Draw("lego same");
    c_FF_pos->cd(2);
    H_FF_XY_light->SetMarkerColor(kRed);  H_FF_XY_light->SetMarkerStyle(kFullCircle);  H_FF_XY_light->SetMarkerSize(0.5); 
    H_FF_XY_heavy->SetMarkerColor(kGreen); H_FF_XY_heavy->SetMarkerStyle(kFullCircle); H_FF_XY_heavy->SetMarkerSize(0.5);
    H_FF_XY_light->Draw("");
    H_FF_XY_heavy->Draw("same");
    c_FF_pos->cd(4);
   	Fission_efficiency->SetAngularOffset(30.);
   	Fission_efficiency->SetEntryRadiusOffset( 4, 0.1);
   	//Fission_efficiency->SetLabelFormat("#splitline{%val (%perc)}{%txt}");
   	Fission_efficiency->SetLabelFormat("#splitline{(%perc)}{%txt}");
   	Fission_efficiency->SetY(.32);
   	Fission_efficiency->Draw("3d t nol");
   	TLegend *pieleg = Fission_efficiency->MakeLegend();
   	pieleg->SetY1(.70); pieleg->SetY2(.90);

/*
    TCanvas *c_FF_angles = new TCanvas("c_FF_angles","Angles of fission fragments",1200,800);
    c_FF_angles->Divide(3,2);
    c_FF_angles->cd(1);
    HTHETA_top_light->SetLineColor(kRed); 	HTHETA_top_heavy->SetLineColor(kGreen);
    HTHETA_top_light->Draw("");				HTHETA_top_heavy->Draw("same");
    c_FF_angles->cd(2);
    HTHETA_bot_light->SetLineColor(kRed); 	HTHETA_bot_heavy->SetLineColor(kGreen);
    HTHETA_bot_light->Draw("");				HTHETA_bot_heavy->Draw("same");
    c_FF_angles->cd(3);
    HTHETA_side_light->SetLineColor(kRed); 	HTHETA_side_heavy->SetLineColor(kGreen);
    HTHETA_side_light->Draw("");			HTHETA_side_heavy->Draw("same");
    c_FF_angles->cd(4);
    HThetaE_top_light->SetMarkerColor(kRed); 		HThetaE_top_heavy->SetMarkerColor(kGreen);
    HThetaE_top_light->SetMarkerStyle(kFullCircle); HThetaE_top_heavy->SetMarkerStyle(kFullCircle);
    HThetaE_top_light->SetMarkerSize(0.5); HThetaE_top_heavy->SetMarkerSize(0.5);
    HThetaE_top_light->Draw("");					HThetaE_top_heavy->Draw("same");
    c_FF_angles->cd(5);
    HThetaE_bot_light->SetMarkerColor(kRed); 		HThetaE_bot_heavy->SetMarkerColor(kGreen);
    HThetaE_bot_light->SetMarkerStyle(kFullCircle);	HThetaE_bot_heavy->SetMarkerStyle(kFullCircle);
    HThetaE_bot_light->SetMarkerSize(0.5);	HThetaE_bot_heavy->SetMarkerSize(0.5);
    HThetaE_bot_light->Draw("");					HThetaE_bot_heavy->Draw("same");
    c_FF_angles->cd(6);
    HThetaE_side_light->SetMarkerColor(kRed); 		HThetaE_side_heavy->SetMarkerColor(kGreen);
    HThetaE_side_light->SetMarkerStyle(kFullCircle);HThetaE_side_heavy->SetMarkerStyle(kFullCircle);
    HThetaE_side_light->SetMarkerSize(0.5);	HThetaE_side_heavy->SetMarkerSize(0.5);
    HThetaE_side_light->Draw("");					HThetaE_side_heavy->Draw("same");


    TCanvas *c_FF_angles1 = new TCanvas("c_FF_angles1","Angles of fission fragments top detector",1200,800);
    c_FF_angles1->cd();
    HThetaE_top_light->SetMarkerColor(kRed); 		HThetaE_top_heavy->SetMarkerColor(kGreen);
    HThetaE_top_light->SetMarkerStyle(kFullCircle); HThetaE_top_heavy->SetMarkerStyle(kFullCircle);
    HThetaE_top_light->SetMarkerSize(0.5); HThetaE_top_heavy->SetMarkerSize(0.5);
    HThetaE_top_light->Draw("colz");					HThetaE_top_heavy->Draw(" col same");
    TCanvas *c_FF_angles2 = new TCanvas("c_FF_angles2","Angles of fission fragments bot detector",1200,800);
    c_FF_angles2->cd();
    HThetaE_bot_light->SetMarkerColor(kRed); 		HThetaE_bot_heavy->SetMarkerColor(kGreen);
    HThetaE_bot_light->SetMarkerStyle(kFullCircle);	HThetaE_bot_heavy->SetMarkerStyle(kFullCircle);
    HThetaE_bot_light->SetMarkerSize(0.5);	HThetaE_bot_heavy->SetMarkerSize(0.5);
    HThetaE_bot_light->Draw("colz");					HThetaE_bot_heavy->Draw("col same");
    TCanvas *c_FF_angles3 = new TCanvas("c_FF_angles3","Angles of fission fragments side detector",1200,800);
    c_FF_angles3->cd();
    HThetaE_side_light->SetMarkerColor(kRed); 		HThetaE_side_heavy->SetMarkerColor(kGreen);
    HThetaE_side_light->SetMarkerStyle(kFullCircle);HThetaE_side_heavy->SetMarkerStyle(kFullCircle);
    HThetaE_side_light->SetMarkerSize(0.5);	HThetaE_side_heavy->SetMarkerSize(0.5);
    HThetaE_side_light->Draw("colz");					HThetaE_side_heavy->Draw("col same");
*/
    TCanvas *c_FF_xy1 = new TCanvas("c_FF_xy1","Position of FF top detector ",1200,800);
    c_FF_xy1->cd();
    H_FF_XY_light_top->Draw("colz");					H_FF_XY_heavy_top->Draw(" col same");
    TCanvas *c_FF_xy2 = new TCanvas("c_FF_xy2","Position of FF bot detector ",1200,800);
    c_FF_xy2->cd();
    H_FF_XY_light_bot->Draw("colz");					H_FF_XY_heavy_bot->Draw(" col same");
    TCanvas *c_FF_xy3 = new TCanvas("c_FF_xy3","Position of FF side detector ",1200,800);
    c_FF_xy3->cd();
    H_FF_XY_light_side->Draw("colz");					H_FF_XY_heavy_side->Draw(" col same");
/*
    TCanvas *c_FF_energyPos1 = new TCanvas("c_FF_energyPos1","Energy per position of fission fragments top detector",1200,800);
    c_FF_energyPos1->cd();
    HyE_top_light->SetMarkerColor(kRed); 		HyE_top_heavy->SetMarkerColor(kGreen);
    HyE_top_light->SetMarkerStyle(kFullCircle);	HyE_top_heavy->SetMarkerStyle(kFullCircle);
    HyE_top_light->SetMarkerSize(0.5);			HyE_top_heavy->SetMarkerSize(0.5);
    HyE_top_light->Draw("colz");				HyE_top_heavy->Draw("col same");
    TCanvas *c_FF_energyPos2 = new TCanvas("c_FF_energyPos2","Energy per position of fission fragments bot detector",1200,800);
    c_FF_energyPos2->cd();
    HyE_bot_light->SetMarkerColor(kRed); 		HyE_bot_heavy->SetMarkerColor(kGreen);
    HyE_bot_light->SetMarkerStyle(kFullCircle);	HyE_bot_heavy->SetMarkerStyle(kFullCircle);
    HyE_bot_light->SetMarkerSize(0.5);			HyE_bot_heavy->SetMarkerSize(0.5);
    HyE_bot_light->Draw("colz");				HyE_bot_heavy->Draw("col same");
    TCanvas *c_FF_energyPos3 = new TCanvas("c_FF_energyPos3","Energy per position of fission fragments side deStector",1200,800);
    c_FF_energyPos3->cd();
    //HxE_side_light->SetMarkerColor(kRed); 		//HxE_side_heavy->SetMarkerColor(kGreen);
    //HxE_side_light->SetMarkerStyle(kFullCircle);//HxE_side_heavy->SetMarkerStyle(kFullCircle);
    //HxE_side_light->SetMarkerSize(0.5);			//HxE_side_heavy->SetMarkerSize(0.5);
    HxE_side_light->Draw("colz");				HxE_side_heavy->Draw("col same");

	//HThetaE_side_heavy->Draw("same");

    TCanvas *c_FF_t = new TCanvas("c_FF_t","Angles of fission fragments",1200,800);
    c_FF_t->cd();
    H_FF_DT_top->SetLineColor(kViolet);
    H_FF_DT_bot->SetLineColor(kCyan);
    H_FF_DT_side->SetLineColor(kOrange);
    H_FF_DT_bot->Draw();
    H_FF_DT_top->Draw("same");
    
    H_FF_DT_side->Draw("same");
	auto legend_velo = new TLegend(0.1,0.75,0.25,0.9);
 	legend_velo->AddEntry(H_FF_DT_top,"Top","f");	
 	legend_velo->AddEntry(H_FF_DT_bot,"Bottom","f");
    legend_velo->AddEntry(H_FF_DT_side,"Side","f");
    legend_velo->Draw();

    TCanvas *c_Ex_pos_HR = new TCanvas("c_Ex_pos_HR","Excitation energy vs position of HR",1200,800);
    c_Ex_pos_HR->Divide(2,1);
    c_Ex_pos_HR->cd(1);
    H_HR_Ex_vs_x->Draw("colz");
    c_Ex_pos_HR->cd(2);
	H_Singles_Eexc_VSTRIP_3->Draw();
*/

    TCanvas *FF_energies = new TCanvas("FF_energies","Energies of fission fragments",1200,800);
    FF_energies->cd();
    H_FF_E_top->SetLineColor(kGreen);
 	H_FF_E_bot->SetLineColor(kBlue);
 	H_FF_E_side->SetLineColor(kRed);
    H_FF_E_top->Draw();
    H_FF_E_bot->Draw("same");
    H_FF_E_side->Draw("same");


  
	
    TFile* File_dp = new TFile("Histograms_dp_fisU239_test.root","RECREATE");
	H_DE_Etot->Write();

	H_FF_XY_light_top->Write();
	H_FF_XY_heavy_top->Write();
	H_FF_XY_light_bot->Write();
	H_FF_XY_heavy_bot->Write();
	H_FF_XY_light_side->Write();
	H_FF_XY_heavy_side->Write();
	c_FF_xy1->Write();
	c_FF_xy2->Write();
	c_FF_xy3->Write();
	Fission_efficiency->Write();
	c_FF_pos->Write();
	
    //for(Int_t j=0;j<6;j++){
    //	H_DE_E[j]->Write();
    //}
//
    //for(Int_t j=0;j<5;j++){
    //	H_E_E[j]->Write();
    //}
//
    //for(Int_t j=0;j<16;j++){
    //	h2_HRmap_Ex_strip[j]->Write();
    //}
    //H_DE_DT->Write();
	//H_Singles_Ekin_VSTRIP_3->Write();
	//H_Singles_Eexc_VSTRIP_3->Write();
	File_dp->Close();
	

	//SW_deut_tree->Write();
	//SW_deut_file->Close();

}



  

     






