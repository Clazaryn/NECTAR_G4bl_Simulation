#include <iostream>
#include <fstream>
#include <stdio.h>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <clocale>

using namespace std;

int GEF_Reader(int Zf, int Af, double Estar, int factor) {

    setlocale(LC_NUMERIC, "C");

    TRandom3* rand = new TRandom3();
    
    //double Estar = (Estar_x10*1.0)/10.0;

    std::ifstream file(Form("../GEF/out/GEFResults_Z%d_A%d_E%.1f_factor_%d.lmd", Zf, Af, Estar, factor));

    if (!file) {
        std::cerr << "File does not exist or could not be opened." << std::endl;
        return 1;  // Exit with error
    }

    if (gSystem->AccessPathName("GEF_tree")) gSystem->mkdir("GEF_tree", kTRUE);

    TFile* outfile = new TFile(Form("./GEF_tree/GEFResults_Z%d_A%d_E%.1f_factor_%d.root", Zf, Af, Estar, factor),"recreate"); //!!!!!!!!!!!!!! 90
    TTree* outputtree = new TTree("GEFtree", "GEF Output in a Root format");

    // Let make the list of the different values that we need 
    int FissionMode;
    Float_t Qvalue; 
    int Z_light; 
    int Z_heavy; 
    int A_light_pre; 
    int A_heavy_pre; 
    int A_light_post; 
    int A_heavy_post; 
    Float_t I1pre; 
    Float_t I2pre; 
    Float_t I1gs; 
    Float_t I2gs; 
    Float_t T_light_pre; 
    Float_t cos_theta_light; 
    Float_t phi_light_deg; 
    Float_t T_heavy_pre; 
    Float_t cos_theta_heavy; 
    Float_t phi_heavy_deg; 
    Float_t Estar_light; 
    Float_t Estar_heavy; 
    int n_light; 
    int n_heavy; 
    Float_t T_light_post; 
    Float_t T_heavy_post; 
    Float_t Estar_fission; 
    Float_t buffer_Float;
    int buffer_int;
    int buffer_int2;

    outputtree->Branch("FissionMode",&FissionMode,"FissionMode/I");
    outputtree->Branch("Qvalue",&Qvalue,"Qvalue/F"); 
    outputtree->Branch("Z_light",&Z_light,"Z_light/I"); 
    outputtree->Branch("Z_heavy",&Z_heavy,"Z_heavy/I"); 
    outputtree->Branch("A_light_pre",&A_light_pre,"A_light_pre/I"); 
    outputtree->Branch("A_heavy_pre",&A_heavy_pre,"A_heavy_pre/I"); 
    outputtree->Branch("A_light_post",&A_light_post,"A_light_post/I"); 
    outputtree->Branch("A_heavy_post",&A_heavy_post,"A_heavy_post/I"); 
    outputtree->Branch("I1pre",&I1pre,"I1pre/F"); 
    outputtree->Branch("I2pre",&I2pre,"I2pre/F"); 
    outputtree->Branch("I1gs",&I1gs,"I1gs/F"); 
    outputtree->Branch("I2gs",&I2gs,"I2gs/F"); 
    outputtree->Branch("T_light_pre",&T_light_pre,"T_light_pre/F"); 
    outputtree->Branch("cos_theta_light",&cos_theta_light,"cos_theta_light/F"); 
    outputtree->Branch("phi_light_deg",&phi_light_deg,"phi_light_deg/F"); 
    outputtree->Branch("T_heavy_pre",&T_heavy_pre,"T_heavy_pre/F"); 
    outputtree->Branch("cos_theta_heavy",&cos_theta_heavy,"cos_theta_heavy/F"); 
    outputtree->Branch("phi_heavy_deg",&phi_heavy_deg,"phi_heavy_deg/F"); 
    outputtree->Branch("Estar_light",&Estar_light,"Estar_light/F"); 
    outputtree->Branch("Estar_heavy",&Estar_heavy,"Estar_heavy/F"); 
    outputtree->Branch("n_light",&n_light,"n_light/I"); 
    outputtree->Branch("n_heavy",&n_heavy,"n_heavy/I"); 
    outputtree->Branch("T_light_post",&T_light_post,"T_light_post/F"); 
    outputtree->Branch("T_heavy_post",&T_heavy_post,"T_heavy_post/F"); 
    outputtree->Branch("Estar_fission",&Estar_fission,"Estar_fission/F"); 


    // NOW WE READ THE FILE

    std::string line;
    int line_num = 0;

    // Step 2: Skip the first 58 lines
    while (line_num < 59 && std::getline(file, line)) {
        line_num++;
    }

    while (std::getline(file, line)) {
    //while (std::getline(file, line) && line_num < 70) {
        // Step 5: Parse the line with sscanf

           //to print the progression
        /*if((line_num)%(5000)==0){
            double j = round(((double)line_num/(double)1000000)*100);     
	    //double j = round(((double)i/(double)10000000)*100);							//!!!!!
            cout << "\r                  LOOP Progression   " << j << " %";
        }*/


        //std::cout << line << std::endl;

        
    

        int result = sscanf(
            line.c_str(),
            "%d %d %d %f %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f",
            &FissionMode, &buffer_int, &buffer_int2, &Qvalue, &Z_light, &Z_heavy, &A_light_pre, &A_heavy_pre, &A_light_post, &A_heavy_post, 
            &I1pre, &I2pre, &I1gs, &I2gs, &T_light_pre, &cos_theta_light, &phi_light_deg, &T_heavy_pre, 
            &cos_theta_heavy, &phi_heavy_deg, &Estar_light, &Estar_heavy, &n_light, &n_heavy, &T_light_post, &T_heavy_post, &Estar_fission
        );

        //if (result != 27) std::cerr << "BAD result=" << result << " line=[" << line << "]\n";

        // Ensure the line matches the expected format (27 values)
        if (result == 27) {
            // Print or process the values (for demonstration, we print them)
            /*std::cout << "Mode: " << FissionMode
                      << " Qvalue: " << Qvalue << " Z1: " << Z_light << " Z2: " << Z_heavy
                      << " A1sci: " << A_light_pre << " A2sci: " << A_heavy_pre << " Elab1pre: " << T_light_pre
                      << " Elab2pre: " << T_heavy_pre << " E@fission: " << Estar_fission << std::endl;*/
            if((line_num)%(5000)==0){
                double j = round(((double)line_num/(double)1000000)*100);     						//!!!!!
                cout << "\r                  LOOP Progression   " << j << " %";
            }

            outputtree->Fill();
            line_num++;
        }

        
        /*else {
            std::cerr << "Error: Line format does not match the expected structure." << std::endl;
        }*/
    }

    outputtree->Write();

    outfile->Close();

}
