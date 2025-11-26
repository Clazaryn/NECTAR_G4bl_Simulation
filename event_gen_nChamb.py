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
import os.path                      
import sys
import datetime
import random
import scipy.constants as const     # this package is for obtaining physical constants
import configparser

from nuclidedatamaster import nuclide_data #this requires the nuclidedatamaster directory
from func_ah import *               # to import file containing functions, eg trigonometric
from deEexc_nChamb import *          # to import file containing deexcitation modes
from nuclear_mass_converter import atomic_to_nuclear_mass  # to convert atomic masses to nuclear masses

# '''''''''' Input parameters ''''''''''''
Nions = int(sys.argv[1])          # number of ions to simulate
HR_type = sys.argv[2]             # heavy recoil type to simulate
excit_En = float(sys.argv[3])     # excitation energy in MeV
ex_label = sys.argv[4]            # excitation energy label

# '''''''''' Reaction parameters ''''''''''''
config = configparser.ConfigParser(inline_comment_prefixes=(';',))
config.read('reac_info_nChamb.txt')

reaction = config['beam_info']['reaction']
A0 = int(config['beam_info']['beam_A'])
qB = int(config['beam_info']['beam_q'])
Z0 = int(config['beam_info']['beam_Z'])
EuA = float(config['beam_info']['beam_EuA'])
sig_pup = float(config['beam_info']['sig_pup'])
emi_x = float(config['beam_info']['emit_x'])
emi_y = float(config['beam_info']['emit_y'])
beta_x = float(config['beam_info']['beta_x'])
beta_y = float(config['beam_info']['beta_y'])
Dispx = float(config['beam_info']['disp_x'])
AT = int(config['target_info']['targ_A'])
ZT = int(config['target_info']['targ_Z'])
Rt = float(config['target_info']['targ_radius'])
x0T = float(config['target_info']['targ_x0'])
AE = int(config['ejectile_info']['ejec_A'])
ZE = int(config['ejectile_info']['ejec_Z'])
AR = int(config['recoil_info']['recoil_A'])
ZR = int(config['recoil_info']['recoil_Z'])
qR = int(config['recoil_info']['recoil_q'])
Sn = float(config['separation_energies']['rec_Sn'])
S2n = float(config['separation_energies']['rec_S2n'])
S3n = float(config['separation_energies']['rec_S3n'])

verbose = False  # boolean for print output

# Starting time to get execution time at the end 
start=datetime.datetime.now()

# '''''''''' Creating and formatting output events files ''''''''''''
# Recoil nucleus & beam
#file_recoil = open("../Event_output/excEn_{}/output_event_generator_recoil_excEn{}.txt".format(ex_label, ex_label), "w")
#file_recoil.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
#file_beam = open("../Event_output/excEn_{}/output_event_generator_beam_excEn{}.txt".format(ex_label, ex_label), "w")
#file_beam.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))

# Gamma emission nucleus
# Ejectile is now created with recoil to ensure pairing in separate files for correlation
if(HR_type == "HRg"):
    file_eject = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HRg_excEn{}_ejectile.txt".format(reaction, ex_label, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_gamma = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HRg_excEn{}_recoil.txt".format(reaction, ex_label, reaction, ex_label), "w")
    file_HR_gamma.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))    

# Neutron emission
elif(HR_type == "HR1n"):
    file_eject = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HR1n_excEn{}_ejectile.txt".format(reaction, ex_label, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HR1n_excEn{}_recoil.txt".format(reaction, ex_label, reaction, ex_label), "w")
    file_HR_neutron.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
elif(HR_type == "HR2n"):
    file_eject = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HR2n_excEn{}_ejectile.txt".format(reaction, ex_label, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron_double = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HR2n_excEn{}_recoil.txt".format(reaction, ex_label, reaction, ex_label),"w")
    file_HR_neutron_double.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))  
elif(HR_type == "HR3n"):
    file_eject = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HR3n_excEn{}_ejectile.txt".format(reaction, ex_label, reaction, ex_label), "w")
    file_eject.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))
    file_HR_neutron_triple = open("../Event_output/{}/excEn_{}/output_event_generator_{}_HR3n_excEn{}_recoil.txt".format(reaction, ex_label, reaction, ex_label),"w")
    file_HR_neutron_triple.write("%s %s %s %s %s %s %s %s %s %s %s %s\n" % ('#x','y','z','px','py','pz','t','PDGid','EventID','TrackID','ParentID','Weight'))  

# _______________________________________________________________________________________________
# REACTION PARAMETERS ----

# '''''''''' Load nuclei data ''''''''''''
amu=931.49410242                        # atomic mass constant from 2025 CODATA
AR_g=AR                                 # heavy residue mass for gamma emission
AR_n=AR-1                               # heavy residue mass for neutron emission
AR_2n=AR-2                              # heavy residue mass for double neutron emission
AR_3n=AR-3                              # heavy residue mass for triple neutron emission

# Masses in MeV (converted from atomic to nuclear masses)
m1 = atomic_to_nuclear_mass(nuclide_data.weight(Z0,A0), Z0)*amu     # Beam 
m2 = atomic_to_nuclear_mass(nuclide_data.weight(ZT,AT), ZT)*amu     # Target
m3 = atomic_to_nuclear_mass(nuclide_data.weight(ZE,AE), ZE)*amu     # Ejectile 
m4 = atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_g), ZR)*amu   # Heavy residue gamma
m5 = atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_n), ZR)*amu   # Heavy residue neutron
m6 = atomic_to_nuclear_mass(nuclide_data.weight(ZR,AR_2n), ZR)*amu  # Heavy residue 2 neutron
m_neutron = const.value(u'neutron mass energy equivalent in MeV')

#E1=EuA*A0
E1=EuA*m1/amu                           # Beam energy in MeV

if(verbose):
    print ("")
    print ("     ___ NUCLEI PARAMETERS ")
    print ("    /")
    print ("-- Number of reaction events : ",Nions)
    print ("-- Beam : A0 =",A0,"; Z0 = ",Z0,"; mass (MeV) = ",m1)
    print ("          Kinetic energy (MeV)/(AMeV) :",E1,"/",EuA)
    print ("-- Target : A1 =",AT,"; Z1 = ",ZT,"; mass (MeV) = ",m2)
    print ("-- Ejectile : A2 =",AE,"; Z2 = ",ZE,"; mass (MeV) = ",m3)
    print ("-- Recoil : A3 =",AR,"; Z3 = ",ZR,"; mass (MeV) = ",m4)
    print ("            Sn (MeV) =",Sn,"; S2n = ",S2n)
    print ("")

# '''''''''' Beam conditions ''''''''''''
sig_x = math.sqrt(beta_x*emi_x)*1000      # rms value of x in mm, it's multiplied by 1000 to have it in mm
sig_y = math.sqrt(beta_y*emi_y)*1000      # rms value of y in mm
sig_ax = math.sqrt(emi_x/beta_x)          # rms value of x' in rad
sig_ay = math.sqrt(emi_y/beta_y)          # rms value of y' in rad

# Global beam/target shift
xshift=0
yshift=0
ztarget=0

if(verbose):
    print ("     ___ BEAM PARAMETERS ")
    print ("    /")
    print("-- Rms position spread (mm) : sigma X =",sig_x)
    print("                              sigma Y = =",sig_y)
    print("-- Rms velocity spread (mm) : sigma aX =",sig_ax)
    print("                              sigma aY =",sig_ay)
    print("-- Rms momentum spread (%) : sigma dP/P =",sig_pup)
    print("-- Emittance (m.rad) : Horizontal = ",emi_x,"; Vertical = ",emi_y)
    print("-- Beta function (m) : Horizontal = ",beta_x,"; Vertical = ",beta_y)
    print("-- Dispersion (m) =",Dispx)

# '''''''''' Variables for plotting results ''''''''''''
plot_x_beam=[]
plot_y_beam=[]
plot_px_beam=[]
plot_py_beam=[]
plot_pz_beam=[]

plot_ejectile_energy=[]
plot_ejectile_theta=[]
plot_ejectile_theta_CM=[]

plot_recoil_energy=[]
plot_recoil_theta=[]
plot_recoil_theta_CM=[]
plot_hr_g_theta=[]
plot_hr_g_energy=[]
plot_hr_n_theta=[]
plot_hr_n_energy=[]
plot_hr_2n_theta=[]
plot_hr_2n_energy=[]
# Control
plot_xslopeP1=[]

# _______________________________________________________________________________________________
# G4BEAMLINE PARTICLE ID CONVERSION ----

# Function to convert Z et A as g4beamline understands it (Z_in and A_in are integer values, 
# Z_out and A_out are values converted for g4beamline) 
    
# '''''''''' PDGid beam  ''''''''''''
# expected 1000922380
Z_beam=int(Z0);
A_beam=int(A0);
PDGid_beam=define_PDGid(Z_beam,A_beam)

# '''''''''' PDGid ejectile ''''''''''''
# expected 1000010010  for (d,p), 1000010020 for (d,d)
Z_eject=int(ZE);
A_eject=int(AE);
PDGid_eject=define_PDGid(Z_eject,A_eject)

# '''''''''' PDGid compound nucleus  ''''''''''''
# expected 1000922390  for (d,p), 1000922380 for (d,d)
Z_recoil=int(ZR);
A_recoil=int(AR);
PDGid_recoil=define_PDGid(Z_recoil,A_recoil)

# '''''''''' PDGid HR gamma emission  ''''''''''''
# expected 1000922390  for (d,p), 1000922380 for (d,d)
Z_HR_g=int(ZR);
A_HR_g=int(AR_g);
PDGid_HR_g=define_PDGid(Z_HR_g,A_HR_g)

# '''''''''' PDGid HR neutron emission  ''''''''''''
# expected 1000922380  for (d,p), 1000922370 for (d,d)
Z_HR_n=int(ZR);
A_HR_n=int(AR_n);
PDGid_HR_n=define_PDGid(Z_HR_n,A_HR_n)

# '''''''''' PDGid HR double neutron emission  ''''''''''''
# expected 1000922370  for (d,p), 1000922360 for (d,d)
Z_HR_2n=int(ZR);
A_HR_2n=int(AR_2n);
PDGid_HR_2n=define_PDGid(Z_HR_2n,A_HR_2n)

# '''''''''' PDGid HR triple neutron emission  ''''''''''''
# expected 1000922370  for (d,p), 1000922360 for (d,d)
Z_HR_3n=int(ZR);
A_HR_3n=int(AR_3n);
PDGid_HR_3n=define_PDGid(Z_HR_3n,A_HR_3n)

# _______________________________________________________________________________________________
# MAIN LOOP ----

if(verbose):
    print ("")
    print ("******* GENERATE EVENTS: run 1 *******")
    print ("")

Ncounter=0
Nreaction=1                     # current reaction event number that fulfills the reaction condition specified line 351
Nbeam=1                         # current ion beam event number (Nbeam>=Nreaction as not each beam ion will interact with the target)
Ekres=0                         # kinetic energy heavy residue
phi=0
N_2n_events=0                   # number of 2 neutron events created
N_3n_events=0                   # number of 3 neutron events created

# /!\ Ncounter is updated/defined at bottom of while loop /!\
while Ncounter < Nions+1:

    # Generation of recoil determined from HR_type
    if(HR_type == "HRg"):
        gamma_in=1
        neutron_in=0
    elif HR_type in ("HR1n", "HR2n", "HR3n"):
        gamma_in=0
        neutron_in=1
    else:
        print("HR_type not recognised: recoil could not be assigned")

    # '''''''''' Excitation energy ''''''''''''
    Eex = excit_En      # excitation energy defined by input
    #if (gamma_in==1):        
    #    Eex=random.uniform(0, Max_eexc)
    #elif(neutron_in==1):
    #    Eex=random.uniform(Sn, Max_eexc)
    Eex_min_2sol=0.4   

    # '''''''''' Interaction coordinates  ''''''''''''
  
    # Beam event : random (x,y) point within emittance and random momentum spread value
    xb = np.random.normal(0,sig_x)
    yb = np.random.normal(0,sig_y)
    xslope = np.random.normal(0,sig_ax)
    yslope = np.random.normal(0,sig_ay)
    delta_pupb = np.random.normal(0,sig_pup)

    # Beam event coordinates at z=0
    x=xb+xshift+Dispx*1000*delta_pupb
    y=yb+yshift
    z=0
    plot_x_beam.append(x)
    plot_y_beam.append(y)

    # Beam momentum
    P1=momentum_rel(E1,m1) 
    px_beam=xslope*P1
    py_beam=yslope*P1
    pz_beam=delta_pupb*P1+P1
    plot_px_beam.append(px_beam)
    plot_py_beam.append(py_beam)
    plot_pz_beam.append(pz_beam)
    
    # Checks and file updates
    # print("Ekin beam : ",E1, " Momentum beam : ",P1, " Nreaction : ", Nreaction, " px_beam : ", px_beam, " py_beam : ", py_beam, " delta_pz_beam : ", delta_pupb*P1)   
    #file_beam.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,px_beam,py_beam,pz_beam,0,int(PDGid_beam),Nbeam,1,0,1))        
    Nbeam=Nbeam+1
    
    # Random z within target diameter
    z_target=np.random.uniform(-Rt,Rt)    
    z=z_target

    # '''''''''' Two body collision and deexcitation chanels  ''''''''''''
    if((math.pow(x-x0T,2)+math.pow(z,2))<math.pow(Rt,2)): # if beam event is in target volume
     
        # Calling two body collision function from deexcitation file
        E3,E4,pX_4,pY_4,pZ_4,pX_3,pY_3,pZ_3,theta_4,theta_3,theta_4_CM,theta_3_CM,phi_4,thetaOK=two_body_col(E1,m1,m2,m3,m4,Eex,AE,Eex_min_2sol)
        #print("theta HR :",theta_4," -  theta_4 : ",theta_4*180/math.pi)
        # Exit loop if ejectile angle is over the max kinematic angle
        if(thetaOK==0):
            #print("Eex :",Eex)
            continue

        # Checks and file updates        
        #print("pX_eject : ",pX_3," pY_eject : ",pY_3," pZ_eject : ",pZ_3)        
        #print("ptot eject : ",math.sqrt(math.pow(pX_3,2)+math.pow(pY_3,2)+math.pow(pZ_3,2)))      
        #file_recoil.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-pX_4,pY_4,pZ_4,0,int(PDGid_recoil),Nreaction,1,0,1))
        
        # '''''''''' Entering here if user chose gamma emission  ''''''''''''
        if(gamma_in==1):
            
            # Calling gamma channel from deexcitation file
            Ek_residue_g = dexec_gamma(m4,qR,E4,Eex,-1)            
            m_residue_g=m4
            #print("theta_4",theta_4*180/math.pi)             
            P_residue_g=momentum_rel(Ek_residue_g,m_residue_g)                       
            px_residue_g=P_residue_g*math.sin(theta_4)*math.cos(phi_4)+xslope*P1 #xslope initial horizontal slope
            py_residue_g=P_residue_g*math.sin(theta_4)*math.sin(phi_4)+yslope*P1 #yslope initial vertical slope
            pz_residue_g=P_residue_g*math.cos(theta_4)+delta_pupb*P1           
            #print("px_residue_g :",px_residue_g,"py_residue_n :",py_residue_g,"pz_residue_g :",pz_residue_g)
            #print("theta_4 :",theta_4,"phi_4 :",phi_4,"theta_3",theta_3,"E3",E3)
            
            P_residue_g_slope=Mag(px_residue_g,py_residue_g,pz_residue_g)
            # Plotting             
            hr_E_g=math.sqrt(math.pow(m_residue_g,2)+math.pow(P_residue_g_slope,2))-m_residue_g
            #hr_E_g=math.pow(Mag(px_residue_g,py_residue_g,pz_residue_g),2)/(2*m_residue_g)
            hr_phi_g,hr_theta_g=get_theta_phi(px_residue_g,py_residue_g,pz_residue_g)
            #print("hrE :",hrE)
            #print("hrphi :",hrphi,"hrtheta",hrtheta*180/math.pi)
            #if(hrphi<0):
            #    hrphi=hrphi+2.*math.pi

            #print("Eex :",Eex," -  theta_3 : ",theta_3*180/math.pi, "E3 : ", E3)
            #print("Ek_residue_g :",Ek_residue_g," -  theta_4 : ",theta_4*180/math.pi)
            # Checks and file updates 
            #if(theta_4*180/math.pi>0.38):
            if(HR_type == "HRg"):
                file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-pX_3,pY_3,pZ_3,0,int(PDGid_eject),Nreaction,1,0,1))
                file_HR_gamma.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,px_residue_g,py_residue_g,pz_residue_g,0,int(PDGid_HR_g),Nreaction,1,0,1))            
            plot_hr_g_theta.append(hr_theta_g)
            plot_hr_g_energy.append(hr_E_g)

        # '''''''''' Entering here if user chose neutron emission  ''''''''''''                
        if(neutron_in==1):

            # Calling neutron channel from deexcitation file
            Ek_residue_n,m_residue_n,Eex_after_n,vcm=dexec_neutron(m4,qR,AR,E4,Eex,Sn,m_neutron)         
            P_residue_n=momentum_rel(Ek_residue_n,m_residue_n)
            #print("theta_4",theta_4*180/math.pi)
            px_residue_n=P_residue_n*math.sin(theta_4)*math.cos(phi_4)+xslope*P1 #xslope initial horizontal slope
            py_residue_n=P_residue_n*math.sin(theta_4)*math.sin(phi_4)+yslope*P1 #yslope initial vertical slope
            pz_residue_n=P_residue_n*math.cos(theta_4)+delta_pupb*P1
            plot_xslopeP1.append(E4)
            #print(" P_residue_n",px_residue_n)
            # First neutron emission - saving HR energy and angles
            hr_E_n=math.sqrt(math.pow(m_residue_n,2)+math.pow(P_residue_n,2))-m_residue_n
            hr_phi_n,hr_theta_n=get_theta_phi(px_residue_n,py_residue_n,pz_residue_n)

            #if(hr_phi<0):
            #    hr_phi=hr_phi+2.*math.pi

            # Checks and file updates 
            #if(theta_4*180/math.pi>0.38):   
            if(HR_type == "HR1n"):
                file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-pX_3,pY_3,pZ_3,0,int(PDGid_eject),Nreaction,1,0,1))
                file_HR_neutron.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,px_residue_n,py_residue_n,pz_residue_n,0,int(PDGid_HR_n),Nreaction,1,0,1))
            plot_hr_n_theta.append(hr_theta_n)
            plot_hr_n_energy.append(hr_E_n)

            # '''''''''' Entering here if double neutron emission is possible  ''''''''''''
            if(Eex_after_n>S2n):

                # Calling neutron channel from deexcitation file - depends on first neutron emission HR
                N_2n_events=N_2n_events+1
                Ek_residue_2n,m_residue_2n,Eex_after_2n,vcm=dexec_neutron(m5,qR,AR_n,Ek_residue_n,Eex_after_n,S2n,m_neutron)
                P_residue_2n=momentum_rel(Ek_residue_2n,m_residue_2n)
                px_residue_2n=P_residue_2n*math.sin(hr_theta_n)*math.cos(hr_phi_n)+xslope*P1 #xslope initial horizontal slope
                py_residue_2n=P_residue_2n*math.sin(hr_theta_n)*math.sin(hr_phi_n)+yslope*P1 #yslope initial vertical slope
                pz_residue_2n=P_residue_2n*math.cos(hr_theta_n)+delta_pupb*P1

                hr_E_2n=math.pow(Mag(px_residue_2n,py_residue_2n,pz_residue_2n),2)/(2*m_residue_2n)
                hr_phi_2n,hr_theta_2n=get_theta_phi(px_residue_2n,py_residue_2n,pz_residue_2n)
                #if(hrphi1<0):
                #    hrphi1=hrphi1+2.*math.pi
                
                #if(theta_4*180/math.pi>0.38):
                # Checks and file updates 
                if(HR_type == "HR2n"):
                    file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-pX_3,pY_3,pZ_3,0,int(PDGid_eject),Nreaction,1,0,1))
                    file_HR_neutron_double.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,px_residue_2n,py_residue_2n,pz_residue_2n,0,int(PDGid_HR_2n),Nreaction,1,0,1))
                plot_hr_2n_theta.append(hr_theta_2n)
                plot_hr_2n_energy.append(hr_E_2n)   

                if(Eex_after_2n>S3n):
                    #print(Eex_after_2n)
                    # Calling neutron channel from deexcitation file - depends on first neutron emission HR
                    N_3n_events=N_3n_events+1
                    Ek_residue_3n,m_residue_3n,Eex_after_3n,vcm=dexec_neutron(m6,qR,AR_2n,Ek_residue_2n,Eex_after_2n,S3n,m_neutron)
                    P_residue_3n=momentum_rel(Ek_residue_3n,m_residue_3n)
                    px_residue_3n=P_residue_3n*math.sin(hr_theta_2n)*math.cos(hr_phi_2n)+xslope*P1 #xslope initial horizontal slope
                    py_residue_3n=P_residue_3n*math.sin(hr_theta_2n)*math.sin(hr_phi_2n)+yslope*P1 #yslope initial vertical slope
                    pz_residue_3n=P_residue_3n*math.cos(hr_theta_2n)+delta_pupb*P1

                    hr_E_3n=math.pow(Mag(px_residue_3n,py_residue_3n,pz_residue_3n),2)/(2*m_residue_3n)
                    hr_phi_3n,hr_theta_3n=get_theta_phi(px_residue_3n,py_residue_3n,pz_residue_3n)
                    #if(hrphi1<0):
                    #    hrphi1=hrphi1+2.*math.pi
                    
                    # Checks and file updates 
                    if(HR_type == "HR3n"):
                        file_eject.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" % (x,y,z,-pX_3,pY_3,pZ_3,0,int(PDGid_eject),Nreaction,1,0,1))
                        file_HR_neutron_triple.write("%.6f %.6f %.2f %.5f %.5f %.5f %i %i %i %i %i %i\n" %(x,y,z,px_residue_3n,py_residue_3n,pz_residue_3n,0,int(PDGid_HR_3n),Nreaction,1,0,1))
                    #plot_hr_2n_theta.append(hr_theta_2n)
                    #plot_hr_2n_energy.append(hr_E_2n)   

        # '''''''''' Update plotting variables  ''''''''''''
        if (Eex-0.1<Eex<Eex+0.1):
        #if ((Eex>0)):
            plot_recoil_theta.append(theta_4*180/math.pi)
            plot_recoil_theta_CM.append(theta_4_CM)
            plot_recoil_energy.append(E4)
            plot_ejectile_theta.append(theta_3*180/math.pi)
            plot_ejectile_theta_CM.append(theta_3_CM)
            plot_ejectile_energy.append(E3)
        
        # '''''''''' Execution percentage  ''''''''''''
        stride = Nions*10/100
        if(Nreaction%stride==0 and verbose):
            percent=int(100 * Nreaction / Nions)
            print("Events completed : {",percent,"} %") 

        Nreaction=Nreaction+1
        if(Nreaction>1e6): continue

    # define Ncounter based on reaction running
    if(HR_type == "HRg"):  Ncounter = Nreaction
    elif(HR_type == "HR1n"): Ncounter = Nreaction
    elif(HR_type == "HR2n"): Ncounter = N_2n_events
    elif(HR_type == "HR3n"): Ncounter = N_3n_events
    else: Ncounter = Nbeam; print("HR_type not recognised: Ncounter assigned to Nbeam")

# Closing the outputfile
#file_beam.close()
#file_recoil.close()
file_eject.close()
if(HR_type == "HRg"):    file_HR_gamma.close()
elif(HR_type == "HR1n"): file_HR_neutron.close()
elif(HR_type == "HR2n"): file_HR_neutron_double.close()
elif(HR_type == "HR3n"): file_HR_neutron_triple.close()
else: print("HR_type not recognised: files not closed")

#print ("")
#print ("Execution time (h:m:s) = ", datetime.datetime.now()-start)    
#print ("")       
#print ("******* CHOSEN KINEMATIC PLOTS *******")
#print ("")
#print ("==>> Close plot window to continue")
    
# _______________________________________________________________________________________________
# PLOTS ----

#fig_beam, axs_beam = plt.subplots(2, 2,figsize=(12, 10))
#fig_beam.suptitle('Beam position and momentum')
#axs_beam[0, 0].scatter(plot_x_beam,plot_y_beam, s=80, marker='+',c='red')
#axs_beam[0, 0].set_title('Position at target (z=0)')
#axs_beam[0, 0].set(xlabel='x beam (mm)', ylabel='y beam (mm)')
#axs_beam[0, 0].set_xlim([-15.0,15.0])
#axs_beam[0, 0].set_ylim([-15.0,15.0])
#axs_beam[0, 1].hist(plot_px_beam, normed=False, bins=100, edgecolor='black',color='red')
#axs_beam[0, 1].set_title('Momentum px (MeV/c)')
#axs_beam[1, 0].hist(plot_py_beam, normed=False, bins=100, edgecolor='black',color='red')
#axs_beam[1, 0].set_title('Momentum py (MeV/c)')
#axs_beam[1, 1].hist(plot_pz_beam, normed=False, bins=100, edgecolor='black',color='red')
#axs_beam[1, 1].set_title('Momentum pz (MeV/c)')
#
#fig_kin, axs_kin = plt.subplots(2, 2,figsize=(12, 10))
#fig_kin.suptitle('Ejectile and recoil kinematics')
#axs_kin[0, 0].scatter(plot_ejectile_theta,plot_ejectile_energy, s=80, marker='o',c='green')
#axs_kin[0, 0].set_title('Ejectile nucleus energy')
#axs_kin[0, 0].set(xlabel='Theta ejectile (°)', ylabel='Kin. energy ejectile (MeV)')
#axs_kin[0, 1].scatter(plot_recoil_theta,plot_recoil_energy, s=80, marker='o', c='orange')
#axs_kin[0, 1].set_title('Recoil nucleus energy')
#axs_kin[0, 1].set(xlabel='Theta recoil (°)', ylabel='Kin. energy recoil (MeV)')
#axs_kin[0, 1].set_xlim([0.0,0.5])
#axs_kin[1, 0].scatter(plot_ejectile_theta,plot_ejectile_theta_CM, s=80, marker='d',c='green')
#axs_kin[1, 0].set_title('Ejectile nucleus angle')
#axs_kin[1, 0].set(xlabel='Theta ejectile LAB (°)', ylabel='Theta ejectile CM (°)')
#axs_kin[1, 1].scatter(plot_recoil_theta,plot_recoil_theta_CM, s=80, marker='d',c='orange')
#axs_kin[1, 1].set_title('Recoil nucleus angle')
#axs_kin[1, 1].set(xlabel='Theta recoil LAB (°)', ylabel='Theta recoil CM (°)')
#
#plt.show()
