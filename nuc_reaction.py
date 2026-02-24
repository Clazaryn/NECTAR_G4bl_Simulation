#!/usr/bin/env python3.5
# --
# NECTAR EXPERIMENT E028 - SIMULATIONS
# Last modification : C. Berthelot - 28/03/2025
# Adapted from : dexec_2.py Michele Sguazzin 

#  ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.
#
#       Relativistic kinematics 
#
#  ~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.

# _______________________________________________________________________________________________
# LIBRARIES, INPUTS, OUTPUT FILES ----

# '''''''''' Libraries ''''''''''''
from __future__ import print_function
import numpy as np
import math
from scipy.stats import gamma
import sys
sys.path.insert(0, 'UtilityScripts')
import phys_functions as funcs

amu = 931.49410242              # atomic mass constant from 2025 CODATA

max_emission_mode = False       # boolean to switch between max recoil cone and realistic multi-step emission

# _______________________________________________________________________________________________
# TWO BODY COLLISION FUNCTION - E1, E3, etc. are kinetic energies

def two_body_col(E1, m1, m2, m3, m4, Eex, particle, Eex_min_2sol, det_setup='new'):
    thetaOK = False     # boolean to indicate if the ejectile angle is kinematically allowed
    
    # Sampling angle: PoP telescope = 51-69 deg; new setup = 0-90 deg (forward hemisphere)
    det_setup_lower = det_setup.strip().lower() if isinstance(det_setup, str) else 'new'
    if det_setup_lower == 'pop':
        theta3_lab = np.random.uniform(51.0, 69.0) * math.pi / 180      # just sample PoP telescope acceptance
    else:
        theta3_lab = np.random.uniform(0.0, 90.0) * math.pi / 180       # new setup covers full forward hemisphere
    
    # Q value       - $$$$ not used in current build $$$$
    # Q = Q_val_proj_targ_ejec(m1, m2, m3, m4)

    # Relativistic kinematics equations - solving for a defined angle theta3_lab of the ejectile
    Etot = E1 + m1 + m2     # total energy in the lab frame
    m4_new = m4 + Eex       # new mass of the target nucleus after excitation
    P1_lab = funcs.momentum_rel(E1, m1)      # beam momentum in the lab frame

    # Centre of Mass quantities - needed to determine solution space
    Beta_CM_vec = np.array([0.0, 0.0, P1_lab]) / Etot                   # CM beta vector is conserved in collision
    Gamma_CM = 1. / math.sqrt(1 - np.dot(Beta_CM_vec, Beta_CM_vec))     # Center of mass final velocity - gamma
    Etot_CM = math.sqrt((m1 + m2)**2 + 2*E1*m2)                      # total energy in the CM frame (had to hack E1/m1 to get MeV/u)
    E3_CM = (Etot_CM**2 + m3**2 - m4_new**2) / (2. * Etot_CM)           # E3 in the CM frame
    Gamma3_ratio = ((E3_CM/m3) / Gamma_CM)**2                           # ratio of Lorentz factor of E3 and CM frame
    # if Gamma3_ratio >= 1, then only the first kinematic solution is allowed. If Gamma3_ratio < 1, then 2 solutions emerge.

    # '''''''''' Solving kinematic quadratic equation for lab frame solutions ''''''''''''
    # next we solve the quadratic equation grouped as parameters of the ejectile energy E3. 
    # This solution is derived in: https://personalpages.surrey.ac.uk/w.catford/kinematics/relativistic-kinematics-carl-wheldon.pdf
    a = 4. * P1_lab**2 * math.cos(theta3_lab)**2 - 4. * Etot**2
    b = 4. * Etot * (Etot**2 - P1_lab**2  + m3**2 - m4_new**2)
    c = (2. * P1_lab**2 * Etot**2 - 2. * m3**2 * Etot**2 + 2. * m3**2 * P1_lab**2 
        + 2. * m4_new**2 * Etot**2 - 2. * m4_new**2 * P1_lab**2 + 2. * m3**2 * m4_new**2 
        - Etot**4. - P1_lab**4. - m3**4. - m4_new**4. 
        - 4. * m3**2 * P1_lab**2 * math.cos(theta3_lab)**2
    )
    delta = b**2 - 4. * a * c

    # Sorting kinematic solutions for ejectile
    if delta < 0:       # no real solutions, exiting function
        return thetaOK, 0, 0, 0, 0, 0, 0

    if Gamma3_ratio >= 1:       # only the first kinematic solution exists
        if (math.cos(theta3_lab) >= 0):             # minus solution is correct for forward hemisphere
            E3_lab = (-b - math.sqrt(delta)) / (2*a) - m3
        else:                                       # plus solution is correct for backward hemisphere
            E3_lab = (-b + math.sqrt(delta)) / (2*a) - m3
        
        if (E3_lab > 0):        # energies must be positive
            thetaOK = True      # valid solution, continue 
        else:
            return thetaOK, 0, 0, 0, 0, 0, 0        # no valid solution, exiting function

    elif (0 < Gamma3_ratio < 1 and theta3_lab < math.pi/2):    # when Gamma3_ratio is strictly less than 1, two kinematic solutions exist
        E3_1sol = (-b + math.sqrt(delta)) / (2*a) - m3
        E3_2sol = (-b - math.sqrt(delta)) / (2*a) - m3
        if (E3_1sol > 0 and E3_2sol > 0):            # energies must be positive and 2 solutions exist
            thetaOK = True      # valid solution, continue 
        else:
            return thetaOK, 0, 0, 0, 0, 0, 0        # no valid solution, exiting function
        
        #chosen_sol = 1
        chosen_sol = np.random.randint(1, 3)       # if random is desired
        if (chosen_sol == 1):       # if chosen_sol = 1, use the first solution
            E3_lab = E3_1sol
        elif(chosen_sol == 2 and Eex > Eex_min_2sol):   # if chosen_sol = 2, use the second solution if it is above the threshold
            E3_lab = E3_2sol
        else:       # if below threshold, set to 0
            return thetaOK, 0, 0, 0, 0, 0, 0                # below threshold, exiting function

    else:
        if Gamma3_ratio < 0: 
            raise ValueError("Unphysical Gamma3_ratio: Gamma3_ratio = ", Gamma3_ratio, " How did we get here?")
        return thetaOK, 0, 0, 0, 0, 0, 0                # below threshold, exiting function
        
    # '''''''''' defining collision product kinematics ''''''''''''
    P3_lab = funcs.momentum_rel(E3_lab, m3)                     # Ejectile momentum in lab frame
    P3x_lab =  P3_lab * math.sin(theta3_lab)                    # momentum components from conservation
    P3z_lab =  P3_lab * math.cos(theta3_lab)

    P4x_lab = -P3x_lab                                          # conservation defines recoil components
    P4z_lab =  P1_lab - P3z_lab
    P4_lab = math.sqrt(P4x_lab**2 + P4z_lab**2)
    theta4_lab = math.atan2(P4x_lab, P4z_lab)                   # recoil angle (robust, correct quadrant)

    E4_lab = Etot - E3_lab - m3 - m4_new                        # Recoil kinetic energy in lab frame
    P4_lab_check = funcs.momentum_rel(E4_lab, m4_new)           # Recoil momentum in lab frame
    if abs(P4_lab - P4_lab_check) > 1e-6*P4_lab:
        raise ValueError("Recoil momentum mismatch: P4_lab = ", P4_lab, " P4_lab_check = ", P4_lab_check)
    
    # set random phi angles
    if det_setup_lower == 'pop':
        phi3 = math.pi + np.random.uniform(-math.pi/20, math.pi/20) # sample for PoP telescope acceptance - sampled around detector position at phi = 180 degrees
    else:
        phi3 = np.random.uniform(0, math.pi)                    # new setup is only simulated above the beam since its symmetric
        
    phi4 = -phi3                                                # Recoil lab phi angle
    if(phi4 < 0):
        phi4 = phi4 + 2 * math.pi
    
    # Momenta in four-vector format
    FV3_lab = np.concatenate(([E3_lab + m3], funcs.spherical_to_cartesian(P3_lab, theta3_lab, phi3)))       # Ejectile four-vector
    FV4_lab = np.concatenate(([E4_lab + m4_new], funcs.spherical_to_cartesian(P4_lab, theta4_lab, phi4)))   # Recoil four-vector
    # CM Angles for plotting
    FV3_CM = funcs.lorentz_boost(FV3_lab, -Beta_CM_vec, Gamma_CM)
    theta3_CM = math.acos(FV3_CM[3] / np.linalg.norm(FV3_CM[1:4]))
    theta4_CM = math.pi - theta3_CM 

    return thetaOK, FV3_lab, theta3_lab, theta3_CM, FV4_lab, theta4_lab, theta4_CM

# _______________________________________________________________________________________________
# GAMMA DEEXCITATION FUNCTION : A relativistic deexcitation function with options for single-step or cascade gamma emission

def deexec_gamma(A_recoil, FV4, Eex):
    # determine CM frame velocity from recoil before decay
    Beta_CM_vec = FV4[1:4] / FV4[0]                                 # beta vector = momentum / total_energy
    Gamma_CM = 1. / math.sqrt(1 - np.linalg.norm(Beta_CM_vec)**2)   # CM velocity - gamma
    FV4_CM = funcs.lorentz_boost(FV4, -Beta_CM_vec, Gamma_CM)       # boost recoil to the CM frame to start the deexcitation process

    # de-excitation function chooses mode of emission based on max_emission_mode
    if max_emission_mode or Eex < 2.0:              # also single step for low excitation energies
        E_gamma = Eex                               # single step gamma emission - worse case scenario

        # Heavy residue momentum (after gamma emission) in the CM frame
        P_gamma_CM = E_gamma                                                  # HR recoil kick is just gamma momentumn from relativity
        theta_gamma_CM = math.acos(np.random.uniform(-1., 1.))             # CM theta isotropic
        phi_gamma_CM = np.random.uniform(0., 2 * math.pi)                  # CM phi isotropic
        FV4_gamma_CM = np.concatenate(([E_gamma], funcs.spherical_to_cartesian(P_gamma_CM, theta_gamma_CM, phi_gamma_CM)))    # construct gamma four-vector in CM frame
        FV4_CM -= FV4_gamma_CM    # subtract gamma four-vector from recoil four-vector to get the heavy residue four-vector in the CM frame

    else:
        sparam = math.pow(A_recoil, 0.4) / 50        # very rough estimate of slope parameter for the gSF
        def gSF_sample(Emax):                        # simple exponential between 0 and Eex for the gSF - not correct but better than uniform
            return math.log(1 + (math.exp(sparam * Emax) - 1) * np.random.uniform(0.0, 1.0)) / sparam

        Eex_new = Eex
        for i in range(3):                          # assuming the entire cascade is always 3 gammas
            if i == 2:
                E_gamma = Eex_new                   # last gamma always goes to the ground state
            else:
                E_gamma = gSF_sample(Eex_new)       # first 2 gammas are emitted following the gSF
            Eex_new = Eex_new - E_gamma             # update the excitation energy for the next iteration

            # Heavy residue momentum (after gamma emission) in the CM frame
            P_gamma_CM = E_gamma                                               # HR recoil kick is just gamma momentumn from relativity
            theta_gamma_CM = math.acos(np.random.uniform(-1., 1.))             # CM theta isotropic
            phi_gamma_CM = np.random.uniform(0., 2 * math.pi)                  # CM phi isotropic
            FV4_gamma_CM = np.concatenate(([E_gamma], funcs.spherical_to_cartesian(P_gamma_CM, theta_gamma_CM, phi_gamma_CM)))    # construct gamma four-vector in CM frame
            FV4_CM -= FV4_gamma_CM    # subtract gamma four-vector from recoil four-vector to get the heavy residue four-vector in the CM frame

    # return the heavy residue momentum in the LAB frame
    return funcs.lorentz_boost(FV4_CM, Beta_CM_vec, Gamma_CM) 

# _______________________________________________________________________________________________
# NEUTRON DEEXCITATION FUNCTION : A relativistic deexcitation function for single-step or Maxwellian distributed single neutron emission

def deexec_neutron(A_recoil, FV4, Eex, mass_neutron, mult, Sn_CN, mass_HR1n, Sn_1nDght=0, mass_HR2n=0, Sn_2nDght=0, mass_HR3n=0, Sn_3nDght=0, mass_HR4n=0):
    if Eex < Sn_CN:
        raise ValueError("Excitation energy is less than the separation energy: Eex = ", Eex, " Sn_CN = ", Sn_CN)
    if mult < 1 or mult > 4:
        raise ValueError("deexec_neutron: Invalid multiplicity: mult = ", mult)

    Eex_new = Eex
    Sn_array = [Sn_CN, Sn_1nDght, Sn_2nDght, Sn_3nDght]
    mass_array = [mass_HR1n, mass_HR2n, mass_HR3n, mass_HR4n]

    # determine CM frame velocity from recoil before decay
    Beta_CM_vec = FV4[1:4] / FV4[0]         # beta vector = momentum / total_energy
    Gamma_CM = 1. / math.sqrt(1 - np.linalg.norm(Beta_CM_vec)**2)   # CM velocity - gamma
    FV4_CM = funcs.lorentz_boost(FV4, -Beta_CM_vec, Gamma_CM)       # boost recoil to the CM frame to start the deexcitation process

    # assume very basic Maxwellian distribution with a temperature depending on the excitation energy
    def maxwellian_sample(Emax, T):                                    
        if T <= 0 or Emax <= 0:
            raise ValueError("deexec_neutron: Temperature or energy maximum is less than or equal to 0: T = ", T, "Emax = ", Emax)

        if Emax / T < 1e-3:       #if available energy is tiny, sample uniformly
            En = Emax * np.random.uniform(0.0, 0.8) # truncating at 0.8 to avoid negative excitation energy    
            return En
        else:
            u_max = gamma.cdf(Emax, a=2.0, scale=T) # CDF at truncation
            u = np.random.uniform(0.0, u_max)       # Sample uniformly in [0, u_max]
            return gamma.ppf(u, a=2.0, scale=T)     # Invert CDF

    for i in range(mult):
        # de-excitation function chooses mode of emission based on max_emission_mode
        if max_emission_mode:
            if i < mult-1:
                En = 0.01                                               # in worst case scenario, first neutrons have tiny kinetic energy
            elif i == mult-1:
                En = max(0.0, Eex_new - Sn_CN - (mult-1)*0.01)          # single step neutron emission - worst case scenario
            else:
                raise ValueError("Invalid number of iterations: i = ", i, " mult = ", mult)
        else:
            Eexc_available = Eex_new - sum(Sn_array[i:mult])            # available energy for the next neutron
            if Eexc_available <= 0:
                raise RuntimeError("deexec_neutron: Insufficient excitation energy for forced multiplicity")
            # to get max neutron energy, we need to include the recoil kinetic energy of the heavy residue
            # solve
            En_max = (Eexc_available * (Eexc_available + 2.0 * mass_array[i])) / (2.0 * (mass_array[i] + mass_neutron + Eexc_available))  
            T_n = math.sqrt(Eex_new / ((A_recoil - i)/8))               # rough estimate of the neutron temperature in MeV
            En = maxwellian_sample(En_max, T_n)

        # Heavy residue momentum (after neutron emission) in the CM frame
        P_neutron_CM = funcs.momentum_rel(En, mass_neutron)             # neutron momentum in the CM frame
        theta_neutron_CM = math.acos(np.random.uniform(-1., 1.))        # CM theta isotropic
        phi_neutron_CM = np.random.uniform(0., 2 * math.pi)             # CM phi isotropic
        FV4_neutron_CM = np.concatenate(([En + mass_neutron], funcs.spherical_to_cartesian(P_neutron_CM, theta_neutron_CM, phi_neutron_CM)))    # construct neutron four-vector in CM frame
        FV4_CM -= FV4_neutron_CM    # subtract neutron four-vector from recoil four-vector to get the heavy residue four-vector in the CM frame

        E_recoil = math.sqrt(P_neutron_CM**2 + mass_array[i]**2) - mass_array[i]
        Eex_new = max(0.0, Eex_new - Sn_array[i] - En - E_recoil)       # update the excitation energy for the next iteration (recoil kinetic energy included)

    # boost recoil kick back to the LAB frame   
    FV4_new = funcs.lorentz_boost(FV4_CM, Beta_CM_vec, Gamma_CM)    # new four-vector after neutron emission

    return FV4_new, Eex_new