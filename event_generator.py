#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --
# NECTAR EXPERIMENT E028 - SIMULATIONS
# Last modification : G. Leckenby - 8/07/2025
# Adapted from : test_event.py by Ana Henriques, Michele Sguazzin and MatLab event generator by Manfred Grieser 

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
import matplotlib.pyplot as plt 
import math                   
import sys
import time
import os
import scipy.constants as const     # this package is for obtaining physical constants
import configparser

from NuclideDataMaster import nuclide_data #this requires the nuclidedatamaster directory
import nuc_reaction          # to import file containing deexcitation modes
import sys
sys.path.insert(0, 'UtilityScripts')
import phys_functions as funcs

verbose = False          # boolean for print + graph output
writeBeamFile = False   # boolean for if to write Beam + Recoil before deexcitation file

# '''''''''' Input parameters ''''''''''''
# These inputs parameters are passed as command line arguments when running the script
Nions = int(sys.argv[1])          # number of ions to simulate
HR_type = sys.argv[2]             # heavy recoil type to simulate
excit_En = float(sys.argv[3])     # excitation energy in MeV
ex_label = sys.argv[4]            # excitation energy label

# '''''''''' Reaction parameters ''''''''''''
# These reaction parameters are read from the reac_info.txt file
config = configparser.ConfigParser(inline_comment_prefixes=(';',))
config.read('reac_info.txt')

reaction = config['beam_info']['reaction']
A0 = int(config['beam_info']['beam_A'])
qB = int(config['beam_info']['beam_q'])
Z0 = int(config['beam_info']['beam_Z'])
EuA = float(config['beam_info']['beam_EuA'])
sig_pup = float(config['beam_info']['sig_pup'])
emi_x = float(config['beam_info']['emit_x'])
emi_y = float(config['beam_info']['emit_y'])
AT = int(config['target_info']['targ_A'])
ZT = int(config['target_info']['targ_Z'])
Rt = float(config['target_info']['targ_radius'])
x0T = float(config['target_info']['targ_x0'])
AE = int(config['ejectile_info']['ejec_A'])
ZE = int(config['ejectile_info']['ejec_Z'])
AR = int(config['recoil_info']['recoil_A'])
ZR = int(config['recoil_info']['recoil_Z'])
qR = int(config['recoil_info']['recoil_q'])
Sn_CN = float(config['separation_energies']['Sn_CN'])
Sn_1nDght = float(config['separation_energies']['Sn_1nDght'])
Sn_2nDght = float(config['separation_energies']['Sn_2nDght'])
Sn_3nDght = float(config['separation_energies']['Sn_3nDght'])
beta_x = float(config['ring_params']['beta_x'])
beta_y = float(config['ring_params']['beta_y'])
Dispx = float(config['ring_params']['disp_x'])

# Starting time to get execution time at the end 
start = time.perf_counter()

# '''''''''' Creating and formatting output events files ''''''''''''
# Create Event_output directory if it doesn't exist
os.makedirs("../{}_sim/Event_output".format(reaction), exist_ok=True)

# Recoil nucleus & beam
if writeBeamFile:
    file_recoil = open("../{}_sim/Event_output/output_event_generator_recoil_excEn{}.txt".format(reaction, ex_label), "w")
    file_recoil.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_beam = open("../{}_sim/Event_output/output_event_generator_beam_excEn{}.txt".format(reaction, ex_label), "w")
    file_beam.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))

# Gamma emission nucleus
# Ejectile is now created with recoil to ensure pairing in separate files for correlation
if(HR_type == "HRg"):
    file_eject = open("../{}_sim/Event_output/output_event_generator_{}_HRg_excEn{}_ejectile.txt".format(reaction, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_gamma = open("../{}_sim/Event_output/output_event_generator_{}_HRg_excEn{}_recoil.txt".format(reaction, reaction, ex_label), "w")
    file_HR_gamma.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))    

# Neutron emission
elif(HR_type == "HR1n"):
    file_eject = open("../{}_sim/Event_output/output_event_generator_{}_HR1n_excEn{}_ejectile.txt".format(reaction, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron = open("../{}_sim/Event_output/output_event_generator_{}_HR1n_excEn{}_recoil.txt".format(reaction, reaction, ex_label), "w")
    file_HR_neutron.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
elif(HR_type == "HR2n"):
    file_eject = open("../{}_sim/Event_output/output_event_generator_{}_HR2n_excEn{}_ejectile.txt".format(reaction, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron_double = open("../{}_sim/Event_output/output_event_generator_{}_HR2n_excEn{}_recoil.txt".format(reaction, reaction, ex_label),"w")
    file_HR_neutron_double.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
elif(HR_type == "HR3n"):
    file_eject = open("../{}_sim/Event_output/output_event_generator_{}_HR3n_excEn{}_ejectile.txt".format(reaction, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron_triple = open("../{}_sim/Event_output/output_event_generator_{}_HR3n_excEn{}_recoil.txt".format(reaction, reaction, ex_label),"w")
    file_HR_neutron_triple.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
elif(HR_type == "HR4n"):
    file_eject = open("../{}_sim/Event_output/output_event_generator_{}_HR4n_excEn{}_ejectile.txt".format(reaction, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron_quadruple = open("../{}_sim/Event_output/output_event_generator_{}_HR4n_excEn{}_recoil.txt".format(reaction, reaction, ex_label),"w")
    file_HR_neutron_quadruple.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight')) 

# _______________________________________________________________________________________________
# REACTION PARAMETERS ----

# '''''''''' Load nuclei data ''''''''''''
amu = 931.49410242                        # atomic mass constant from 2025 CODATA
AR_g = AR                                 # heavy residue mass for gamma emission
AR_n = AR-1                               # heavy residue mass for neutron emission
AR_2n = AR-2                              # heavy residue mass for double neutron emission
AR_3n = AR-3                              # heavy residue mass for triple neutron emission
AR_4n = AR-4                              # heavy residue mass for quadruple neutron emission

# Masses in MeV (converted from atomic to nuclear masses)
mBeam = funcs.atomic_to_nuclear_mass(nuclide_data.weight(Z0,A0), Z0) * amu     # Beam 
mTarg = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZT,AT), ZT) * amu     # Target
mEjec = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZE,AE), ZE) * amu     # Ejectile 
mHRg = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_g), ZR) * amu    # Heavy residue gamma
mHR1n = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_n), ZR) * amu   # Heavy residue neutron
mHR2n = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_2n), ZR) * amu  # Heavy residue 2 neutron
mHR3n = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_3n), ZR) * amu  # Heavy residue 3 neutron
mHR4n = funcs.atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_4n), ZR) * amu  # Heavy residue 4 neutron
m_neutron = const.value(u'neutron mass energy equivalent in MeV')

E1 = EuA * mBeam / amu                           # Beam energy in MeV

if(verbose):
    print ("")
    print ("     ___ NUCLEI PARAMETERS ")
    print ("    /")
    print ("-- Number of reaction events : ",Nions)
    print ("-- Beam : A0 =", A0, "; Z0 = ", Z0, "; mass (MeV) = ", mBeam)
    print ("          Kinetic energy (MeV)/(AMeV) :", E1, "/", EuA)
    print ("-- Target : A1 =", AT, "; Z1 = ",ZT,"; mass (MeV) = ", mTarg)
    print ("-- Ejectile : A2 =", AE, "; Z2 = ", ZE, "; mass (MeV) = ", mEjec)
    print ("-- Recoil : A3 =", AR, "; Z3 = ", ZR, "; mass (MeV) = ", mHRg)
    print ("            Sn_CN (MeV) =", Sn_CN, "; Sn_1nDght = ", Sn_1nDght, "; Sn_2nDght = ", Sn_2nDght)
    print ("")

# '''''''''' Beam conditions ''''''''''''
sig_y = math.sqrt(beta_y*emi_y) * 1000    # rms value of y in mm
sig_x = math.sqrt(beta_x*emi_x) * 1000    # rms value of x in mm, it's multiplied by 1000 to have it in mm
sig_ax = math.sqrt(emi_x/beta_x)          # rms value of x' in rad
sig_ay = math.sqrt(emi_y/beta_y)          # rms value of y' in rad

# Global beam/target shift
xshift = 0
yshift = 0
ztarget = 0

if(verbose):
    print ("     ___ BEAM PARAMETERS ")
    print ("    /")
    print("-- Rms position spread (mm) : sigma X =", sig_x)
    print("                              sigma Y = =", sig_y)
    print("-- Rms velocity spread (mm) : sigma aX =", sig_ax)
    print("                              sigma aY =", sig_ay)
    print("-- Rms momentum spread (%) : sigma dP/P =", sig_pup)
    print("-- Emittance (m.rad) : Horizontal = ", emi_x, "; Vertical = ", emi_y)
    print("-- Beta function (m) : Horizontal = ", beta_x, "; Vertical = ", beta_y)
    print("-- Dispersion (m) =", Dispx)

# '''''''''' Variables for plotting results ''''''''''''
plot_x_beam = []
plot_y_beam = []
plot_px_beam = []
plot_py_beam = []
plot_pz_beam = []

plot_ejectile_energy = []
plot_ejectile_theta = []
plot_ejectile_theta_CM = []

plot_recoil_energy = []
plot_recoil_theta = []
plot_recoil_theta_CM = []

# _______________________________________________________________________________________________
# G4BEAMLINE PARTICLE ID CONVERSION ----
# Function to convert Z et A as g4beamline understands it - check that isotopes are correctly converted in printed output
# (Z_in and A_in are integer values, Z_out and A_out are values converted for g4beamline) 
    
# '''''''''' PDGid beam  '''''''''''' 
Z_beam = int(Z0)
A_beam = int(A0)
PDGid_beam = funcs.define_PDGid(Z_beam, A_beam)

# '''''''''' PDGid ejectile ''''''''''''
Z_eject = int(ZE)
A_eject = int(AE)
PDGid_eject = funcs.define_PDGid(Z_eject, A_eject)

# '''''''''' PDGid compound nucleus  ''''''''''''
Z_recoil = int(ZR)
A_recoil = int(AR)
PDGid_recoil = funcs.define_PDGid(Z_recoil, A_recoil)

# '''''''''' PDGid HR gamma emission  ''''''''''''
Z_HR_g = int(ZR)
A_HR_g = int(AR_g)
PDGid_HR_g = funcs.define_PDGid(Z_HR_g, A_HR_g)

# '''''''''' PDGid HR neutron emission  ''''''''''''
Z_HR_n = int(ZR)
A_HR_n = int(AR_n)
PDGid_HR_n = funcs.define_PDGid(Z_HR_n, A_HR_n)

# '''''''''' PDGid HR double neutron emission  ''''''''''''
Z_HR_2n = int(ZR)
A_HR_2n = int(AR_2n)
PDGid_HR_2n = funcs.define_PDGid(Z_HR_2n, A_HR_2n)

# '''''''''' PDGid HR triple neutron emission  ''''''''''''
Z_HR_3n = int(ZR)
A_HR_3n = int(AR_3n)
PDGid_HR_3n = funcs.define_PDGid(Z_HR_3n, A_HR_3n)

# '''''''''' PDGid HR quadruple neutron emission  ''''''''''''
Z_HR_4n = int(ZR)
A_HR_4n = int(AR_4n)
PDGid_HR_4n = funcs.define_PDGid(Z_HR_4n, A_HR_4n)

# _______________________________________________________________________________________________
# MAIN LOOP ----

if(verbose):
    print ("")
    print ("******* GENERATE EVENTS: run 1 *******")
    print ("")

Ncounter = 0
Nreaction = 1                     # current reaction event number 
Nbeam = 1                         # current ion beam event number (Nbeam >= Nreaction as not each beam ion will interact with the target)

# /!\ Ncounter is updated/defined at bottom of while loop /!\
while Ncounter < Nions+1:

    # '''''''''' Excitation energy ''''''''''''
    Eex = excit_En      # excitation energy defined by input
    Eex_min_2sol = 0.4   

    # '''''''''' Interaction coordinates  ''''''''''''
  
    # Beam event : random (x,y) point within emittance and random momentum spread value
    xb = np.random.normal(0, sig_x)
    yb = np.random.normal(0, sig_y)
    xslope = np.random.normal(0, sig_ax)
    yslope = np.random.normal(0, sig_ay)
    delta_pupb = np.random.normal(0, sig_pup)

    # Beam event coordinates at z=0
    x = xb + xshift + Dispx * 1000 * delta_pupb
    y = yb + yshift
    z = 0
    plot_x_beam.append(x)
    plot_y_beam.append(y)

    # calculate beam properties for plotting and file writing
    P1 = funcs.momentum_rel(E1, mBeam)
    plot_px_beam.append(P1 * xslope)
    plot_py_beam.append(P1 * yslope)
    plot_pz_beam.append(P1 * (1 + delta_pupb))
    
    # Checks and file updates
    # print("Ekin beam : ",E1, " Momentum beam : ",P1, " Nreaction : ", Nreaction, " px_beam : ", px_beam, " py_beam : ", py_beam, " delta_pz_beam : ", delta_pupb*P1)   
    if writeBeamFile: file_beam.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,P1*xslope,P1*yslope,P1*(1+delta_pupb),0,int(PDGid_beam),Nbeam,1,0,1))        
    Nbeam = Nbeam + 1
    
    # Random z within target diameter
    z_target = np.random.uniform(-Rt, Rt)    
    z = z_target

    # '''''''''' Two body collision and deexcitation chanels are implemented in nuc_reaction.py  ''''''''''''
    # in 2-body collision: 1 corresponds to beam, 2 to target, 3 to ejectile, 4 to heavy residue
    if((math.pow(x - x0T,2) + math.pow(z,2)) < math.pow(Rt,2)): # if beam event is in target volume
     
        # Calling two body collision function from nuc_reaction.py
        thetaOK, FV3, theta3, theta3_CM, FV4, theta4, theta4_CM = nuc_reaction.two_body_col(E1, mBeam, mTarg, mEjec, mHRg, Eex, AE, Eex_min_2sol)
        # print("theta HR :",theta4," -  theta4 : ",theta4*180/math.pi)
        
        if not thetaOK:     # skip to next event if ejectile angle is not valid (not kinematically allowed)
            # print("Eex :",Eex)
            continue

        # adds emmittance and momentum spread to the recoil four-vector
        p4_copy = FV4[1:4].copy()      
        p4_copy[0] += p4_copy[2] * xslope                               # angular emittance
        p4_copy[1] += p4_copy[2] * yslope        
        p4_copy[2] *= (1.0 + delta_pupb)                                # momentum spread
        FV4[1:4] = p4_copy                                              # update the recoil four-vector
        FV4[0] = math.sqrt((mHRg+Eex)**2 + np.dot(p4_copy, p4_copy))    # puts four-vector back on mass shell after update (includes excitation energy)

        # Checks and file updates        
        # print("ptot eject : ",math.sqrt(math.pow(pX_3,2)+math.pow(pY_3,2)+math.pow(pZ_3,2)))      
        # print("pX_eject : ",pX_3," pY_eject : ",pY_3," pZ_eject : ",pZ_3)        
        if writeBeamFile: file_recoil.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV4[1],FV4[2],FV4[3],0,int(PDGid_recoil),Nreaction,1,0,1))
        
        # '''''''''' Gamma emission recoil function ''''''''''''
        if(HR_type == "HRg"):            
            FV4_deexc = nuc_reaction.deexec_gamma(AR, FV4, Eex)     # execute gamma decay
            mHR_final = mHRg  # final heavy residue mass after deexcitation

            # write to output files
            file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV3[1],FV3[2],FV3[3],0,int(PDGid_eject),Nreaction,1,0,1))
            file_HR_gamma.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV4_deexc[1],FV4_deexc[2],FV4_deexc[3],0,int(PDGid_HR_g),Nreaction,1,0,1))

        # '''''''''' 1 neutron emission recoil function ''''''''''''
        elif(HR_type == "HR1n"):         
            FV4_deexc, Eex_resid = nuc_reaction.deexec_neutron(AR, FV4, Eex, m_neutron, 1, Sn_CN, mHR1n)     # execute neutron decay
            if Eex_resid > 0.01:
                FV4_deexc = nuc_reaction.deexec_gamma(AR_n, FV4_deexc, Eex_resid)     # execute an additional gamma decay
            mHR_final = mHR1n  # final heavy residue mass after deexcitation

            # write to output files
            file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV3[1],FV3[2],FV3[3],0,int(PDGid_eject),Nreaction,1,0,1))
            file_HR_neutron.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,-FV4_deexc[1],FV4_deexc[2],FV4_deexc[3],0,int(PDGid_HR_n),Nreaction,1,0,1))

        # '''''''''' 2 neutron emission recoil function ''''''''''''
        elif(HR_type == "HR2n"):
            FV4_deexc, Eex_resid = nuc_reaction.deexec_neutron(AR, FV4, Eex, m_neutron, 2, Sn_CN, mHR1n, Sn_1nDght, mHR2n)     # execute neutron decay
            if Eex_resid > 0.01:
                FV4_deexc = nuc_reaction.deexec_gamma(AR_2n, FV4_deexc, Eex_resid)     # execute an additional gamma decay
            mHR_final = mHR2n  # final heavy residue mass after deexcitation

            # write to output files
            file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV3[1],FV3[2],FV3[3],0,int(PDGid_eject),Nreaction,1,0,1))
            file_HR_neutron_double.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,-FV4_deexc[1],FV4_deexc[2],FV4_deexc[3],0,int(PDGid_HR_2n),Nreaction,1,0,1))
               
        # '''''''''' 3 neutron emission recoil function ''''''''''''
        elif(HR_type == "HR3n"):
            FV4_deexc, Eex_resid = nuc_reaction.deexec_neutron(AR, FV4, Eex, m_neutron, 3, Sn_CN, mHR1n, Sn_1nDght, mHR2n, Sn_2nDght, mHR3n)     # execute neutron decay
            if Eex_resid > 0.01:
                FV4_deexc = nuc_reaction.deexec_gamma(AR_3n, FV4_deexc, Eex_resid)     # execute an additional gamma decay
            mHR_final = mHR3n  # final heavy residue mass after deexcitation
               
            # write to output files
            file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV3[1],FV3[2],FV3[3],0,int(PDGid_eject),Nreaction,1,0,1))
            file_HR_neutron_triple.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,-FV4_deexc[1],FV4_deexc[2],FV4_deexc[3],0,int(PDGid_HR_3n),Nreaction,1,0,1))
        
        # '''''''''' 4 neutron emission recoil function ''''''''''''
        elif(HR_type == "HR4n"):
            FV4_deexc, Eex_resid = nuc_reaction.deexec_neutron(AR, FV4, Eex, m_neutron, 4, Sn_CN, mHR1n, Sn_1nDght, mHR2n, Sn_2nDght, mHR3n, Sn_3nDght, mHR4n)     # execute neutron decay
            if Eex_resid > 0.01:
                FV4_deexc = nuc_reaction.deexec_gamma(AR_4n, FV4_deexc, Eex_resid)     # execute an additional gamma decay
            mHR_final = mHR4n  # final heavy residue mass after deexcitation
            # write to output files
            file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-FV3[1],FV3[2],FV3[3],0,int(PDGid_eject),Nreaction,1,0,1))
            file_HR_neutron_quadruple.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,-FV4_deexc[1],FV4_deexc[2],FV4_deexc[3],0,int(PDGid_HR_4n),Nreaction,1,0,1))
        
        else:
            raise ValueError("event_generator: HR_type not recognised: ", HR_type)

        # '''''''''' Update plotting variables  ''''''''''''
        plot_recoil_theta.append(theta4 * 180 / math.pi)
        plot_recoil_theta_CM.append(theta4_CM * 180 / math.pi)
        plot_recoil_energy.append(FV4_deexc[0] - mHR_final)  # kinetic energy of heavy residue after deexcitation
        plot_ejectile_theta.append(theta3 * 180 / math.pi)
        plot_ejectile_theta_CM.append(theta3_CM * 180 / math.pi)
        plot_ejectile_energy.append(FV3[0] - mEjec)  # kinetic energy of ejectile
        
        # '''''''''' Execution percentage  ''''''''''''
        stride = int(Nions * 10 / 100)
        if Ncounter % stride == 0 and verbose:
            percent=int(100 * Ncounter / Nions)
            print("Events completed : {",percent,"} %") 

        Ncounter += 1
        Nreaction += 1

# Closing the outputfile
if writeBeamFile: 
    file_beam.close()
    file_recoil.close()
file_eject.close()
if(HR_type == "HRg"):    file_HR_gamma.close()
elif(HR_type == "HR1n"): file_HR_neutron.close()
elif(HR_type == "HR2n"): file_HR_neutron_double.close()
elif(HR_type == "HR3n"): file_HR_neutron_triple.close()
elif(HR_type == "HR4n"): file_HR_neutron_quadruple.close()
else: print("HR_type not recognised: files not closed")

if verbose:
    print ("")
    print ("Execution time (h:m:s) = ", time.perf_counter()-start)    
    print ("")       
    print ("******* CHOSEN KINEMATIC PLOTS *******")
    print ("")
    print ("==>> Close plot window to continue")
    
    # _______________________________________________________________________________________________
    # PLOTS ----
    fig_beam, axs_beam = plt.subplots(2, 2,figsize=(12, 10))
    fig_beam.suptitle('Beam position and momentum (no target coincidence)')
    axs_beam[0, 0].scatter(plot_x_beam,plot_y_beam, s=10, marker='+',c='red')
    axs_beam[0, 0].set_title('Position at target (z=0)')
    axs_beam[0, 0].set(xlabel='x beam (mm)', ylabel='y beam (mm)')
    axs_beam[0, 0].set_xlim([-15.0,15.0])
    axs_beam[0, 0].set_ylim([-15.0,15.0])
    axs_beam[0, 1].hist(plot_px_beam, normed=False, bins=100, edgecolor='black',color='red')
    axs_beam[0, 1].set_title('Momentum px (MeV/c)')
    axs_beam[1, 0].hist(plot_py_beam, normed=False, bins=100, edgecolor='black',color='red')
    axs_beam[1, 0].set_title('Momentum py (MeV/c)')
    axs_beam[1, 1].hist(plot_pz_beam, normed=False, bins=100, edgecolor='black',color='red')
    axs_beam[1, 1].set_title('Momentum pz (MeV/c)')

    fig_kin, axs_kin = plt.subplots(2, 2,figsize=(12, 10))
    fig_kin.suptitle('Ejectile and recoil kinematics')
    axs_kin[0, 0].scatter(plot_ejectile_theta,plot_ejectile_energy, s=10, marker='o',c='green')
    axs_kin[0, 0].set_title('Ejectile nucleus energy')
    axs_kin[0, 0].set(xlabel='Theta ejectile (°)', ylabel='Kin. energy ejectile (MeV)')
    axs_kin[0, 1].scatter(plot_recoil_theta,plot_recoil_energy, s=10, marker='o', c='orange')
    axs_kin[0, 1].set_title('Recoil nucleus energy')
    axs_kin[0, 1].set(xlabel='Theta recoil (°)', ylabel='Kin. energy recoil (MeV)')
    axs_kin[0, 1].set_xlim([0.0,0.5])
    axs_kin[1, 0].scatter(plot_ejectile_theta,plot_ejectile_theta_CM, s=10, marker='d',c='green')
    axs_kin[1, 0].set_title('Ejectile nucleus angle')
    axs_kin[1, 0].set(xlabel='Theta ejectile LAB (°)', ylabel='Theta ejectile CM (°)')
    axs_kin[1, 1].scatter(plot_recoil_theta,plot_recoil_theta_CM, s=10, marker='d',c='orange')
    axs_kin[1, 1].set_title('Recoil nucleus angle')
    axs_kin[1, 1].set(xlabel='Theta recoil LAB (°)', ylabel='Theta recoil CM (°)')

    plt.show()
