// --
// NECTAR EXPERIMENT E028 - SIMULATIONS
// Last modification : C. Berthelot - 13/12/2023
// Adapted from : Fission1.C created by Michele Sguazzin

//  ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
//
//       				Fission fragments 
//
//  ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

// _________________________________________________________________________________________________
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
#include <TRandom3.h>

void Fission1(int Nevents_recoil, double mass_recoil){

	// _______________________________________________________________________________________________
	// INPUT/OUTPUT FILES, HISTOGRAMS, VARIABLE DECLARATION ----

	// ========= Input files ==========
  char FissionFragments_CM[]="../GEF_output/239U_6to20mev.lmd"; 										// File .lmd input from GEF : Fission framents in the CM frame 
  ifstream FF_CM_stream;
  // NOTE: Update reaction name in path as needed (e.g., "../238U_dp_sim/Event_output/...")
  char RecoilNuclei_Lab[]="../Event_output/output_event_generator_recoil_nofields.txt";	// File .txt input from kinematics : recoil nuclei in the lab frame
 	ifstream Recoil_stream;

	// ========= Output files =========
  // NOTE: Update reaction name in path as needed (e.g., "../238U_dp_sim/Event_output/...")
  char FissionFragments_Lab[]="../Event_output/FissionFragments_Lab.txt";  			// File .txt output : fission fragments in the lab frame
  ofstream FF_Lab_stream;

  // ========= Histograms ===========

  // Angular distribution of fragments   
  TH1D *HThetaLab_heavy = new TH1D("ThetaLab_heavy","Lab angle of heavy fragments; #theta_{Fis,lab} (#circ); Counts",80,0,25);
  TH1D *HThetaLab_light = new TH1D("ThetaLab_light","Lab angle of light fragments; #theta_{Fis,lab} (#circ); Counts",80,0,25);
  TH1D *HThetaCM_heavy = new TH1D("ThetaCM_heavy","CM angle of heavy fragments; #theta^{CM} (deg)",60,0,180);
  TH1D *HThetaCM_light = new TH1D("ThetaCM_light","CM angle of light fragments; #theta^{CM} (deg)",60,0,180);

  // Mass distribution of fragments 
  TH1D *HMass_heavy = new TH1D("HMass_heavy","Heavy fission fragments mass A; Mass A",120,60,180);
  TH1D *HMass_light = new TH1D("HMass_light","Light fission fragments mass A; Mass A",120,60,180);

  // Energy and velocity of fragments 
  TH1D *HveloCM_light = new TH1D("HveloCM_light","Velocity distribution of fission fragments (CM); Velocity (cm/ns)",100,0,2);
  TH1D *HveloCM_heavy = new TH1D("HveloCM_heavy","Velocity distribution of fission fragments (CM); Velocity (cm/ns)",100,0,2);  
  TH1D *Hvelo_light = new TH1D("Hvelo_light","Velocity distribution of fission fragments (Lab); Velocity (cm/ns)",100,3,7.5);
  TH1D *Hvelo_heavy = new TH1D("Hvelo_heavy","Velocity distribution of fission fragments (Lab); Velocity (cm/ns)",100,3,7.5);
	TH1D *HEk_light = new TH1D("HEk_light","Kinetic energy distribution of fission fragments (Lab); Kinetic energy of fission fragments (MeV); Counts",50,0,3500);
  TH1D *HEk_heavy = new TH1D("HEk_heavy","Kinetic energy distribution of fission fragments (Lab); Kinetic energy of fission fragments (MeV); Counts",50,0,3500);
  TH1D *HEk_light_AMeV = new TH1D("HEk_light_AMeV","Kinetic energy distribution of fission fragments (Lab); Kinetic energy of fission fragments (MeV/u); Counts",50,0,30);
  TH1D *HEk_heavy_AMeV = new TH1D("HEk_heavy_AMeV","Kinetic energy distribution of fission fragments (Lab); Kinetic energy of fission fragments (MeV/u); Counts",50,0,30);
 	TH1D *Hvelo_recoil = new TH1D("Hvelo_recoil","Recoil nucleus velocity lab frame; Velocity (cm/ns)",100,5.1,5.2);
 	TH1D *Henergy_recoil = new TH1D("Henergy_recoil","Recoil nucleus kinetic energy lab frame; Ek (MeV)",500,2800,3400);
 	
 	// Energy and velocity vs angle 
 	TH2D *HThetaE_light = new TH2D("HThetaE_light","Fission fragments energy per angle;#theta (deg);Ek (MeV)",120,0,20,4000,0,4000);
  TH2D *HThetaE_heavy = new TH2D("HThetaE_heavy","Fission fragments energy per angle;#theta (deg);Ek (MeV)",120,0,20,4000,0,4000);
 	TH2D *HThetaVelo_light = new TH2D("HThetaVelo_light","Fission fragments velocity per angle;#theta (deg);v (cm/ns)",120,0,20,100,3.5,7.5);
 	TH2D *HThetaVelo_heavy = new TH2D("HThetaVelo_heavy","Fission fragments velocity per angle;#theta (deg);v (cm/ns)",120,0,20,100,3.5,7.5); 	
  
  // Control 
 	TH2D *HThetaCMVlab_light = new TH2D("HThetaCMVlab_light","Fission fragments lab velocity vs CM angle;#theta^{CM} (deg);Velocity (cm/ns)",60,0,180,100,3,7.5); 	
 	TH2D *HThetaCMVlab_heavy = new TH2D("HThetaCMVlab_heavy","Fission fragments lab velocity vs CM angle;#theta^{CM} (deg);Velocity (cm/ns)",60,0,180,100,3,7.5); 
  TH2D *HThetaCMLab_heavy = new TH2D("ThetaCMLab_heavy","Lab angle vs CM angle of fission fragments; #theta^{CM} (deg); #theta^{lab}(deg)",60,0,180,50,0,25);
  TH2D *HThetaCMLab_light = new TH2D("ThetaCMLab_light","Lab angle vs CM angle of fission fragments; #theta^{CM} (deg); #theta^{lab}(deg)",60,0,180,50,0,25);
  
  // ========= Variable declaration =========

  // Recoil input file fields
	double x_in, y_in, z_in,Px_in,Py_in,Pz_in,t_in,PDGid_in,EventID_in,TrackID_in,ParentID_in,Weight_in;

	// Fission fragments input file fields
	int GEFsize=100000;
  double Z_light[GEFsize], Z_heavy[GEFsize], A_light[GEFsize], A_heavy[GEFsize], EkCM_light[GEFsize], EkCM_heavy[GEFsize];

  // Intermediate calculation variables
  double P_recoil, E_recoil, M_recoil, gamma_recoil, beta_x_recoil,beta_y_recoil,beta_z_recoil,v_recoil;				// Recoil nucleus : center of mass
  double M_light,PCM_light,PxCM_light,PyCM_light,PzCM_light,ECM_light,thetaCM_light,vCM_light;									// Light fission fragment
  double P_light,Px_light,Py_light,Pz_light,phi_light,thetaLab_light,Ek_light,v_light; 
  double M_heavy,PCM_heavy,PxCM_heavy,PyCM_heavy,PzCM_heavy,ECM_heavy,thetaCM_heavy,vCM_heavy;									// Heavy fission fragment
  double P_heavy,Px_heavy,Py_heavy,Pz_heavy,phi_heavy,thetaLab_heavy,Ek_heavy, v_heavy; 

  // Pulse height defect variables
  double resistivity=3000; 			// Micron detector resistivity = 3 to 10 kOhm.cm
  double mobility=1481;					// Constant electron mobility
  double biasVoltage=300;				// BB29 bias voltage 
  double chargeCentroid = 0 ;				

  // Fission fragments output file fields
  int Nevents_FF=(int)(Nevents_recoil*2);
  
 	double *outX = new double[Nevents_FF];	double *outY = new double[Nevents_FF];
 	double *outZ = new double[Nevents_FF];	double *outPX = new double[Nevents_FF];
  double *outPY = new double[Nevents_FF];	double *outPZ = new double[Nevents_FF];
 	double *outT = new double[Nevents_FF];
 	int *outPDGID = new int[Nevents_FF];		int *outEventID = new int[Nevents_FF];
 	int *outTrackID = new int[Nevents_FF];	int *outParentID = new int[Nevents_FF];
 	int *outWeight = new int[Nevents_FF];		int *outNZ = new int[Nevents_FF];
 	int *outNA = new int[Nevents_FF];

 	// This type of array declaration limits the size to 40000 elements 
	//double outX[Nevents_FF],outY[Nevents_FF],outZ[Nevents_FF],outPX[Nevents_FF],outPY[Nevents_FF],outPZ[Nevents_FF],outT[Nevents_FF];
  //int outPDGID[Nevents_FF],outEventID[Nevents_FF],outTrackID[Nevents_FF],outParentID[Nevents_FF],outWeight[Nevents_FF],outNZ[Nevents_FF],outNA[Nevents_FF];
 
  // Counters and randomizators 
  int n_FF_in=0, n_recoil=0, n_FF_out=0;
  TRandom3 *random = new TRandom3();

  // ========= Anisotropy functions =========
  // Implementation of the anisotropy effect on the angular distribution of fission fragments
  // It allows to choose between 6 values of anisotropy factor: [0.0,0.2,0.6,1.0,1.5,2.0]
  // By assigning a corresponding value to the 'ani' variable:  [ 0 , 1 , 2 , 3 , 4 , 5 ]

  TF1 *f1 = new TF1("f1","1+0.0*TMath::Power(x,2)",-1,1);
  TF1 *f2 = new TF1("f2","1+0.2*TMath::Power(x,2)",-1,1);
  TF1 *f3 = new TF1("f3","1+0.6*TMath::Power(x,2)",-1,1);
  TF1 *f4 = new TF1("f4","1+1.0*TMath::Power(x,2)",-1,1);
  TF1 *f5 = new TF1("f5","1+1.5*TMath::Power(x,2)",-1,1);
  TF1 *f6 = new TF1("f6","1+2.0*TMath::Power(x,2)",-1,1);
  int ani = 0;	// Choosing ani=0 leads to anisotropy factor of 0.0

	// _______________________________________________________________________________________________
	// FISSION EVENTS SIMULATION ----
       
  // ========= Fission events in CM frame =========
  // This loop gets the fission events from the GEF file and stores the fields into arrrays

  FF_CM_stream.open(FissionFragments_CM);   
  while(true){
		// Breaks if the file is wrong (mainly if first ascii lines with fields description haven't been removed)
    if (!FF_CM_stream.good()) break;		
    // Filling the FF arrays with GEF file
    FF_CM_stream>>Z_light[n_FF_in]>>Z_heavy[n_FF_in]>>A_light[n_FF_in]>>A_heavy[n_FF_in]>>EkCM_light[n_FF_in]>>EkCM_heavy[n_FF_in];
    // Calculating center of mass velocity of fragments
    vCM_light=sqrt((2*EkCM_light[n_FF_in]*1.6*pow(10,-13))/(A_light[n_FF_in]*1.66054*pow(10,-27)))*pow(10,-7);
    vCM_heavy=sqrt((2*EkCM_heavy[n_FF_in]*1.6*pow(10,-13))/(A_heavy[n_FF_in]*1.66054*pow(10,-27)))*pow(10,-7);
    // Filling histograms
    HMass_heavy->Fill(A_heavy[n_FF_in]);	HMass_light->Fill(A_light[n_FF_in]);    
    HveloCM_light->Fill(vCM_light);				HveloCM_heavy->Fill(vCM_heavy);

    n_FF_in++;    
  }
  FF_CM_stream.close();

  // ========= Recoil events in LAB frame =========
  // This loop gets the recoil events from the kinematics and for each recoil nucleus, 
  // two fission events in the lab frame are associated (heavy and light fragments) and stored in the output file
     
  Recoil_stream.open(RecoilNuclei_Lab);     	
  while(true){
  	// Breaks if the file is wrong (mainly if first ascii lines with fields description haven't been removed)
		if (!Recoil_stream.good()) break;
		// Filling recoil variables with file events one by one
		Recoil_stream>>x_in>>y_in>>z_in>>Px_in>>Py_in>>Pz_in>>t_in>>PDGid_in>>EventID_in>>TrackID_in>>ParentID_in>>Weight_in; 
	
		// ========= Recoil nucleus kinematics =========
		// Relativistic kinematics of recoil nucleus (Momentum, energy, gamma and beta velocities)
		P_recoil = sqrt(pow(Px_in,2)+pow(Py_in,2)+pow(Pz_in,2));
		E_recoil = sqrt(pow(P_recoil,2)+pow(mass_recoil,2));
		gamma_recoil = E_recoil/mass_recoil;
		beta_x_recoil = Px_in/E_recoil;	beta_y_recoil = Py_in/E_recoil;	beta_z_recoil = Pz_in/E_recoil;	
		v_recoil = sqrt((2*(E_recoil-mass_recoil)*1.6*pow(10,-13))/((mass_recoil/931.5)*1.66054*pow(10,-27)))*pow(10,-7);
		
		// ========= Fission fragment CM angle sampling =========
		// Samples theta and phi for light and heavy fragments depending on anisotropy value applied
		if(ani==0){	thetaCM_light = acos(f1->GetRandom(-1,1)); }
		if(ani==1){	thetaCM_light = acos(f2->GetRandom(-1,1)); }
		if(ani==2){	thetaCM_light = acos(f3->GetRandom(-1,1)); }
		if(ani==3){	thetaCM_light = acos(f4->GetRandom(-1,1)); }
		if(ani==4){	thetaCM_light = acos(f5->GetRandom(-1,1)); }
		if(ani==5){	thetaCM_light = acos(f6->GetRandom(-1,1)); }   
 
		phi_light = gRandom->Uniform(0,TMath::Pi()*2);
		phi_heavy = phi_light+TMath::Pi();      
		thetaCM_heavy = TMath::Pi()-thetaCM_light;

		if(phi_heavy>TMath::Pi()*2){
			 phi_heavy=phi_heavy-TMath::Pi()*2;
		}
	
		// ========= Fission fragment kinematics =========
		// Calculates relativistic kinematics for each fragment in the CM frame and 
		// adds the relativistic movement of the CM frame (recoil nucleus)

		// Masses in MeV -- !! Not exact masses from nuclear data !! 
		M_light = A_light[n_recoil]*931.5; M_heavy = A_heavy[n_recoil]*931.5; 
				
		// Light fragment kinematics
		// Center of mass frame (CM momentums (x,y,z) and energy)
		PCM_light = sqrt(pow(EkCM_light[n_recoil]+M_light,2)-pow(M_light,2));
		ECM_light = EkCM_light[n_recoil]+M_light;
		PxCM_light = PCM_light*sin(thetaCM_light)*cos(phi_light);
		PyCM_light = PCM_light*sin(thetaCM_light)*sin(phi_light);
		PzCM_light = PCM_light*cos(thetaCM_light);
		// Lab frame (lab momentum (x,y,z), angle, energy and velocity)
		Px_light = gamma_recoil*(PxCM_light+beta_x_recoil*ECM_light);
		Py_light = gamma_recoil*(PyCM_light+beta_y_recoil*ECM_light);
		Pz_light = gamma_recoil*(PzCM_light+beta_z_recoil*ECM_light);
		P_light = sqrt(pow(Px_light,2)+pow(Py_light,2)+pow(Pz_light,2));
		thetaLab_light = acos(Pz_light/P_light);
		Ek_light = sqrt(pow(P_light,2)+pow(M_light,2))-M_light;
		v_light = sqrt((2*Ek_light*1.6*pow(10,-13))/(M_light/931.5*1.66054*pow(10,-27)))*pow(10,-7);

		// Pulse height defect



		// Filling output arrays for light fragments	
		if(n_FF_out<Nevents_FF){
			outX[n_FF_out]=x_in;outY[n_FF_out]=y_in; outZ[n_FF_out]=z_in; 
			outPX[n_FF_out]=Px_light; outPY[n_FF_out]=Py_light; outPZ[n_FF_out]=Pz_light; 
			outT[n_FF_out]=0; outNZ[n_FF_out]=Z_light[n_recoil] ;outNA[n_FF_out]=A_light[n_recoil]; 
			outEventID[n_FF_out]=n_recoil+1; outTrackID[n_FF_out]=1; outParentID[n_FF_out]=0; outWeight[n_FF_out]=1;

			n_FF_out++;
		}
		
		// Heavy fragment kinematics
		// Center of mass frame (CM momentums (x,y,z) and energy)
		PCM_heavy = sqrt(pow(EkCM_heavy[n_recoil]+M_heavy,2)-pow(M_heavy,2));
		ECM_heavy = EkCM_heavy[n_recoil]+M_heavy;
		PxCM_heavy = PCM_heavy*sin(thetaCM_heavy)*cos(phi_heavy);
		PyCM_heavy = PCM_heavy*sin(thetaCM_heavy)*sin(phi_heavy); 
		PzCM_heavy = PCM_heavy*cos(thetaCM_heavy);	
		// Lab frame (lab momentum (x,y,z), angle, energy and velocity)
		Px_heavy = gamma_recoil*(PxCM_heavy+beta_x_recoil*ECM_heavy);
		Py_heavy = gamma_recoil*(PyCM_heavy+beta_y_recoil*ECM_heavy);
		Pz_heavy = gamma_recoil*(PzCM_heavy+beta_z_recoil*ECM_heavy);	
		P_heavy = sqrt(pow(Px_heavy,2)+pow(Py_heavy,2)+pow(Pz_heavy,2));
		thetaLab_heavy = acos(Pz_heavy/P_heavy);
		Ek_heavy = sqrt(pow(P_heavy,2)+pow(M_heavy,2))-M_heavy;
		v_heavy = sqrt((2*Ek_heavy*1.6*pow(10,-13))/(M_heavy/931.5*1.66054*pow(10,-27)))*pow(10,-7);

		// Filling output arrays for heavy fragments
	  if(n_FF_out<Nevents_FF){
			outX[n_FF_out]=x_in;outY[n_FF_out]=y_in; outZ[n_FF_out]=z_in; 
			outPX[n_FF_out]=Px_heavy; outPY[n_FF_out]=Py_heavy; outPZ[n_FF_out]=Pz_heavy; 
			outT[n_FF_out]=0; outNZ[n_FF_out]=Z_heavy[n_recoil]; outNA[n_FF_out]=A_heavy[n_recoil]; 
			outEventID[n_FF_out]=n_recoil+1; outTrackID[n_FF_out]=1; outParentID[n_FF_out]=0; outWeight[n_FF_out]=1;
   		
   		n_FF_out++;
		}
		
		// ========= Filling histograms =========	
		HThetaLab_light->Fill(thetaLab_light*180/TMath::Pi());	HThetaLab_heavy->Fill(thetaLab_heavy*180/TMath::Pi());
    HThetaCM_light->Fill(thetaCM_light*180/TMath::Pi());		HThetaCM_heavy->Fill(thetaCM_heavy*180/TMath::Pi());   
    HThetaCMLab_light->Fill(thetaCM_light*180/TMath::Pi(),thetaLab_light*180/TMath::Pi());	HThetaCMLab_heavy->Fill(thetaCM_heavy*180/TMath::Pi(),thetaLab_heavy*180/TMath::Pi());    
    Hvelo_light->Fill(v_light);	Hvelo_heavy->Fill(v_heavy);    	
    Hvelo_recoil->Fill(v_recoil);	Henergy_recoil->Fill(E_recoil-mass_recoil);
    HThetaE_light->Fill(thetaLab_light*180/TMath::Pi(),Ek_light);	HThetaE_heavy->Fill(thetaLab_heavy*180/TMath::Pi(),Ek_heavy);    
    HEk_light->Fill(Ek_light);	HEk_heavy->Fill(Ek_heavy);
    HEk_light_AMeV->Fill(Ek_light/(M_light/931.5));	HEk_heavy_AMeV->Fill(Ek_heavy/(M_heavy/931.5));      
    HThetaVelo_light->Fill(thetaLab_light*180/TMath::Pi(),v_light);	HThetaVelo_heavy->Fill(thetaLab_heavy*180/TMath::Pi(),v_heavy);    
    HThetaCMVlab_light->Fill(thetaCM_light*180/TMath::Pi(),v_light);	HThetaCMVlab_heavy->Fill(thetaCM_heavy*180/TMath::Pi(),v_heavy);
		
		n_recoil++;
	}  
	Recoil_stream.close();

  // ========= Filling output file =========
  FF_Lab_stream.open(FissionFragments_Lab);

  for(int i=0;i<Nevents_FF;i++){

    if(outNA[i]<100){
      FF_Lab_stream<<outX[i]<<"  "<<outY[i]<<"  "<<outZ[i]<<"  "<<outPX[i]<<"  "<<outPY[i]<<"  "<<outPZ[i]<<"  "<<outT[i]<<"  1000"<<outNZ[i]<<"0"<<outNA[i]<<"0  "<<outEventID[i]<<"  "<<outTrackID[i]<<"  "<<outParentID[i]<<"  "<<outWeight[i]<<"\n";     
    }

    if(outNA[i]>=100){
      FF_Lab_stream<<outX[i]<<"  "<<outY[i]<<"  "<<outZ[i]<<"  "<<outPX[i]<<"  "<<outPY[i]<<"  "<<outPZ[i]<<"  "<<outT[i]<<"  1000"<<outNZ[i]<<outNA[i]<<"0  "<<outEventID[i]<<"  "<<outTrackID[i]<<"  "<<outParentID[i]<<"  "<<outWeight[i]<<"\n"; 
    }  
  }
  FF_Lab_stream.close();

	// _______________________________________________________________________________________________
	// CANVAS DISPLAY ----

  // ========= Angular distibutions =========
  // in center of mass and lab frames 
  auto legend_anisotropy = new TLegend(0.4,0.10,0.6,0.4);
  legend_anisotropy->SetHeader("Anisotropy factor","C");
  legend_anisotropy->AddEntry(f1,"0.0","f");	legend_anisotropy->AddEntry(f2,"0.2","f");
  legend_anisotropy->AddEntry(f3,"0.6","f");	legend_anisotropy->AddEntry(f4,"1.0","f");
  legend_anisotropy->AddEntry(f5,"1.5","f");	legend_anisotropy->AddEntry(f6,"2.0","f");
  auto legend_thetaCM = new TLegend(0.1,0.75,0.25,0.9);
  legend_thetaCM->AddEntry(HThetaCM_heavy,"Heavy fragments","f");	
  legend_thetaCM->AddEntry(HThetaCM_light,"Light fragments","f");
  auto legend_thetalab = new TLegend(0.1,0.75,0.25,0.9);
  legend_thetalab->AddEntry(HThetaLab_heavy,"Heavy fragments","f");	
  legend_thetalab->AddEntry(HThetaLab_light,"Light fragments","f");
  auto legend_tEkinlab = new TLegend(0.1,0.75,0.25,0.9);
  legend_tEkinlab->AddEntry(HEk_heavy,"Heavy fragments","f");	
  legend_tEkinlab->AddEntry(HEk_light,"Light fragments","f");
  auto legend_tEkinlab_AMeV = new TLegend(0.1,0.75,0.25,0.9);
  legend_tEkinlab_AMeV->AddEntry(HEk_heavy_AMeV,"Heavy fragments","f");	
  legend_tEkinlab_AMeV->AddEntry(HEk_light_AMeV,"Light fragments","f");

  TCanvas *ctheta = new TCanvas("ctheta","Fragment angular and mass distributions",1000,800);
  ctheta->Divide(2,2);
  ctheta->cd(1);    
	HThetaCM_heavy->SetLineColor(kGreen);	HThetaCM_light->SetLineColor(kRed);
  HThetaCM_light->SetTitle("Angular distribution of fission fragments (CM)");
  HThetaCM_light->Draw(); HThetaCM_heavy->Draw("same");  
  legend_thetaCM->Draw();    
  ctheta->cd(2);    
  HThetaLab_heavy->SetLineColor(kGreen); HThetaLab_light->SetLineColor(kRed);
  HThetaLab_heavy->SetTitle("Angular distribution of fission fragments (lab)");
  HThetaLab_heavy->Draw(); HThetaLab_light->Draw("same");
  legend_thetaCM->Draw();   
  ctheta->cd(3);
  f1->SetTitle("Sampling Cos(#theta) for different anisotropy factors");
  f1->GetXaxis()->SetTitle("Cos(#theta)");
  f1->SetLineColor(kOrange+7);f1->Draw();
  f2->SetLineColor(kPink+7);f2->Draw("same");
  f3->SetLineColor(kMagenta-2);f3->Draw("same");
  f4->SetLineColor(kAzure+8);f4->Draw("same");
  f5->SetLineColor(kTeal-8);f5->Draw("same");
  f6->SetLineColor(kYellow-5);f6->Draw("same");
	legend_anisotropy->Draw();
	ctheta->cd(4);
  HMass_heavy->SetLineColor(kGreen); 	HMass_light->SetLineColor(kRed);     
  HMass_light->Draw();  							HMass_heavy->Draw("same");
  legend_thetaCM->Draw();

  // ========= Velocity/Energy distibutions ========= 
  // Recoil nucleus in lab frame and FF in center of mass and lab frames
	auto legend_velo = new TLegend(0.1,0.75,0.25,0.9);
  legend_velo->AddEntry(HEk_heavy,"Heavy fragments","f");	
  legend_velo->AddEntry(HEk_light,"Light fragments","f");

  TCanvas *cvelo = new TCanvas("cvelo","Fragment velocity and energy distributions",1000,800);
  cvelo->Divide(2,2);
  cvelo->cd(1);    
  HEk_light->SetLineColor(kRed); 				HEk_heavy->SetLineColor(kGreen); 
  HEk_light->Draw();										HEk_heavy->Draw("same"); 
  legend_velo->Draw();  
  cvelo->cd(2);
  HThetaE_light->SetMarkerColor(kRed);	HThetaE_heavy->SetMarkerColor(kGreen);
  HThetaE_light->Draw();								HThetaE_heavy->Draw("same");
  legend_velo->Draw();
  cvelo->cd(3);
  Hvelo_light->SetLineColor(kRed);			Hvelo_heavy->SetLineColor(kGreen);  
  Hvelo_heavy->Draw();									Hvelo_light->Draw("same");
  legend_velo->Draw();
  cvelo->cd(4);   
  HThetaVelo_light->SetMarkerColor(kRed);	HThetaVelo_heavy->SetMarkerColor(kGreen);
  HThetaVelo_light->Draw();								HThetaVelo_heavy->Draw("same");
  legend_velo->Draw();
							
  // ========= Others ========= 
  // Theta vs energy of fragment, theta vs velocity of fragments 
  TCanvas *ccontrol = new TCanvas("ccontrol","Recoil nucleus and fission fragment CM velocity ",1000,800);
  ccontrol->Divide(2,2);
  ccontrol->cd(1);
  Hvelo_recoil->Draw();  
  ccontrol->cd(2);
  HveloCM_light->SetLineColor(kRed); 	HveloCM_heavy->SetLineColor(kGreen); 
  HveloCM_light->Draw();	HveloCM_heavy->Draw("same"); 
  ccontrol->cd(3);
  HThetaCMLab_light->SetMarkerColor(kRed);HThetaCMLab_light->SetMarkerStyle(kFullCircle);HThetaCMLab_light->SetMarkerSize(0.5);
  HThetaCMLab_heavy->SetMarkerColor(kGreen);HThetaCMLab_heavy->SetMarkerStyle(kFullCircle);HThetaCMLab_heavy->SetMarkerSize(0.5);
  HThetaCMLab_light->Draw("colz");  
  HThetaCMLab_heavy->Draw("col same"); 
  ccontrol->cd(4);
  HThetaCMVlab_light->SetMarkerColor(kRed);HThetaCMVlab_light->SetMarkerStyle(kFullCircle);HThetaCMVlab_light->SetMarkerSize(0.5);
  HThetaCMVlab_heavy->SetMarkerColor(kGreen);HThetaCMVlab_heavy->SetMarkerStyle(kFullCircle);HThetaCMVlab_heavy->SetMarkerSize(0.5);
  HThetaCMVlab_light->Draw("colz");
  HThetaCMVlab_heavy->Draw("col same");

  TCanvas *cangle = new TCanvas("cangle","cangle",1000,800);
  cangle->cd();
  HThetaLab_heavy->SetLineColor(kGreen+2); HThetaLab_light->SetLineColor(kRed);
  HThetaLab_heavy->SetTitle("Angular distribution of fission fragments (lab)");
  HThetaLab_light->Draw(); HThetaLab_heavy->Draw("same");
  legend_thetalab->Draw(); 

  TCanvas *cTkin = new TCanvas("cTkin","cTkin",1000,800);
  cTkin->cd();
  HEk_heavy->SetLineColor(kGreen+2); HEk_light->SetLineColor(kRed);
  HEk_heavy->SetTitle("Kinetic energies of fission fragments (lab)");
  HEk_light->Draw(); HEk_heavy->Draw("same");
  legend_tEkinlab->Draw();     

  TCanvas *cTkin_AMeV = new TCanvas("cTkin_AMeV","cTkin_AMeV",1000,800);
  cTkin_AMeV->cd();
  HEk_heavy_AMeV->SetLineColor(kGreen+2); HEk_light_AMeV->SetLineColor(kRed);
  HEk_heavy_AMeV->SetTitle("Kinetic energies of fission fragments (lab)");
  HEk_heavy_AMeV->Draw(); HEk_light_AMeV->Draw("same");
  legend_tEkinlab_AMeV->Draw();  

}


  

     






