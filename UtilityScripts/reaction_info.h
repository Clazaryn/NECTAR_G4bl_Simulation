#ifndef REACTION_INFO_H
#define REACTION_INFO_H

#include "ini_parser.h"
#include <iostream>
#include <vector>

// Structure to hold reaction information
struct ReactionInfo {
    // Beam info
    std::string reaction;
    int beam_A, beam_q, beam_Z;
    double beam_EuA, sig_pup, emit_x, emit_y, beta_x, beta_y, disp_x;
    
    // Target info
    int targ_A, targ_Z;
    double targ_radius, targ_x0;
    
    // Ejectile info
    int ejec_A, ejec_Z;
    
    // Recoil info
    int recoil_A, recoil_Z, recoil_q;
    std::vector<double> recoil_excEns;
    
    // Separation energies
    double rec_Sn, rec_S2n, rec_S3n;
    
    // Reaction masses
    double mass_beam, mass_target, mass_recoil, mass_eject;
    
    // Ring parameters
    std::string det_setup;  // Detector setup: "new" or "PoP"
    double dip_radius, quad_Leff, dip_Brho, dip_Bfld;
    double quad_Q1, quad_Q2, quad_Q3, quad_Q4, quad_Q5;
    
    // Load all values from INI file
    bool loadFromIni(const std::string& filename) {
        try {
            IniParser parser;
            if (!parser.loadFile(filename)) {
                std::cerr << "Error: Could not load INI file " << filename << std::endl;
                return false;
            }
            
            // Beam info
            reaction = parser.getString("beam_info", "reaction");
            beam_A = parser.getInt("beam_info", "beam_A");
            beam_q = parser.getInt("beam_info", "beam_q");
            beam_Z = parser.getInt("beam_info", "beam_Z");
            beam_EuA = parser.getDouble("beam_info", "beam_EuA");
            sig_pup = parser.getDouble("beam_info", "sig_pup");
            emit_x = parser.getDouble("beam_info", "emit_x");
            emit_y = parser.getDouble("beam_info", "emit_y");
            beta_x = parser.getDouble("beam_info", "beta_x");
            beta_y = parser.getDouble("beam_info", "beta_y");
            disp_x = parser.getDouble("beam_info", "disp_x");
            
            // Target info
            targ_A = parser.getInt("target_info", "targ_A");
            targ_Z = parser.getInt("target_info", "targ_Z");
            targ_radius = parser.getDouble("target_info", "targ_radius");
            targ_x0 = parser.getDouble("target_info", "targ_x0");
            
            // Ejectile info
            ejec_A = parser.getInt("ejectile_info", "ejec_A");
            ejec_Z = parser.getInt("ejectile_info", "ejec_Z");
            
            // Recoil info
            recoil_A = parser.getInt("recoil_info", "recoil_A");
            recoil_Z = parser.getInt("recoil_info", "recoil_Z");
            recoil_q = parser.getInt("recoil_info", "recoil_q");
            recoil_excEns = parser.getDoubleVector("recoil_info", "recoil_excEns");
            
            // Separation energies
            rec_Sn = parser.getDouble("separation_energies", "rec_Sn");
            rec_S2n = parser.getDouble("separation_energies", "rec_S2n");
            rec_S3n = parser.getDouble("separation_energies", "rec_S3n");
            
            // Reaction masses
            mass_beam = parser.getDouble("reaction_masses", "mass_beam");
            mass_target = parser.getDouble("reaction_masses", "mass_target");
            mass_recoil = parser.getDouble("reaction_masses", "mass_recoil");
            mass_eject = parser.getDouble("reaction_masses", "mass_eject");
            
            // Ring parameters
            det_setup = parser.getString("ring_params", "det_setup");
            dip_radius = parser.getDouble("ring_params", "dip_radius");
            quad_Leff = parser.getDouble("ring_params", "quad_Leff");
            dip_Brho = parser.getDouble("ring_params", "dip_Brho");
            dip_Bfld = parser.getDouble("ring_params", "dip_Bfld");
            quad_Q1 = parser.getDouble("ring_params", "quad_Q1");
            quad_Q2 = parser.getDouble("ring_params", "quad_Q2");
            quad_Q3 = parser.getDouble("ring_params", "quad_Q3");
            quad_Q4 = parser.getDouble("ring_params", "quad_Q4");
            quad_Q5 = parser.getDouble("ring_params", "quad_Q5");
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading INI file: " << e.what() << std::endl;
            return false;
        }
    }
    
    // Print all loaded values for debugging
    void printInfo() {
        std::cout << "=== Reaction Information ===" << std::endl;
        std::cout << "Reaction: " << reaction << std::endl;
        std::cout << "Beam: A=" << beam_A << ", Z=" << beam_Z << ", q=" << beam_q << ", E/A=" << beam_EuA << " MeV/u" << std::endl;
        std::cout << "Target: A=" << targ_A << ", Z=" << targ_Z << ", radius=" << targ_radius << " mm" << std::endl;
        std::cout << "Ejectile: A=" << ejec_A << ", Z=" << ejec_Z << std::endl;
        std::cout << "Recoil: A=" << recoil_A << ", Z=" << recoil_Z << ", q=" << recoil_q << std::endl;
        std::cout << "Masses: beam=" << mass_beam << ", target=" << mass_target << ", recoil=" << mass_recoil << ", ejectile=" << mass_eject << " MeV" << std::endl;
        std::cout << "Detector setup: " << det_setup << std::endl;
        std::cout << "Excitation energies: ";
        for (size_t i = 0; i < recoil_excEns.size(); ++i) {
            std::cout << recoil_excEns[i];
            if (i < recoil_excEns.size() - 1) std::cout << ", ";
        }
        std::cout << " MeV" << std::endl;
    }
};

#endif // REACTION_INFO_H 