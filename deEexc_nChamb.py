#!/usr/bin/env python
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
import os.path
import random
import scipy.constants as const                 # this package is for obtaining physical constants

from nuclidedatamaster import nuclide_data      # this requires the nuclidedatamaster directory
from func_ah import *                           # to import a file containing functions
from nuclear_mass_converter import atomic_to_nuclear_mass  # to convert atomic masses to nuclear masses

amu=931.49410242                                # atomic mass constant from 2025 CODATA

# _______________________________________________________________________________________________
# TWO BODY COLLISION FUNCTION 

def two_body_col(E1,m1,m2,m3,m4,Eex,particle,Eex_min_2sol):

    E3_threshold=0.05

    # Sampling angle in the whole telescope area (taking into account target radius angular resolution)
    theta3 = np.arccos(1 - 2 * np.random.uniform(0.0, 0.5))
    #theta3=np.random.uniform(35.0,90.0)*math.pi/180;
    #theta3=60*math.pi/180;
    
    # Q value 
    Q=Q_val_proj_targ_ejec(m1,m2,m3,m4)  
    #print("Excitation energy Eex (MeV) : ",Eex,"Q value : ",Q)
    #print("Theta ejectile -> degree : ",theta3*180/math.pi,"radian : ",theta3,"cosinus : ",math.cos(theta3))

    # Relativistic kinematics equations
    Etot = E1 + m1 + m2
    m4_new=m4+Eex
    P1 = math.sqrt(math.pow(E1,2)+2.*E1*m1)
    a=4.*math.pow(P1,2)*math.pow(math.cos(theta3),2.)-4.*math.pow(Etot,2.)
    b=4.*math.pow(Etot,3.)-4.*math.pow(P1,2.)*Etot+4.*math.pow(m3,2.)*Etot-4.*math.pow(m4_new,2.)*Etot
    c=2.*math.pow(P1,2.)*math.pow(Etot,2.)-2.*math.pow(m3,2.)*math.pow(Etot,2.)+2.*math.pow(m3,2.)*math.pow(P1,2.)+2.*math.pow(m4_new,2.)*math.pow(Etot,2.)-2.*math.pow(m4_new,2.)*math.pow(P1,2.)+2.*math.pow(m3,2.)*math.pow(m4_new,2.)-math.pow(Etot,4.)-math.pow(P1,4.)-math.pow(m4_new,4.)-math.pow(m3,4.)-4.*math.pow(m3,2.)*math.pow(P1,2.)*math.pow(math.cos(theta3),2.)
    delta = math.pow(b,2.)-4.*a*c
    #print("Etot : ",Etot,"Beam momentum P1 : ",P1,"Beam energy E1 : ",E1)
    #print("a",a,"b",b,"c",c,"delta",delta);

    # Calculating kinematic solutions for ejectile
    # Quadratic equation
    if(delta>=0):   # Delta < 0
        E3_m=(-b-math.sqrt(delta))/(2*a)-m3
        E3_p=(-b+math.sqrt(delta))/(2*a)-m3
        if(E3_m>0 and E3_p>0):
            n_sol=2
            thetaOK=1
            #print("okay")       
        else:
            n_sol=0
            E3=0
    elif(delta==0): # Delta = 0
        thetaOK=1
        n_sol=2
        E3_m=-b/(2*a)-m3
        E3_p=-b/(2*a)-m3
    else:           # Delta > 0
        n_sol=0
        E3=0
        thetaOK=0

    # In case of 2 solutions, we take a random one 
    # (and checks the energy threshold for 2nd solution)
    #random_solution=random.randint(1,2)
    random_solution=1
    if(n_sol==2):   # 2 solutions or 1 solution
        if(E3_p>=E3_m):
            E3_1sol=E3_p
            E3_2sol=E3_m
        else:
            E3_1sol=E3_m
            E3_2sol=E3_p
        if (random_solution==1):
            E3=E3_1sol
            #print("E3",E3)
        elif(random_solution==2 and Eex>Eex_min_2sol):
            E3=E3_2sol
            #print("E3",E3)
        else:
            E3=0
            n_sol=0

    if(n_sol==0): # exit loop
        #print("No solution")
        thetaOK=0
        E3,E4,pX_4,pY_4,pZ_4,pX_3,pY_3,pZ_3,theta4,theta3,theta4_CM,theta3_CM,phi4_lab=0,0,0,0,0,0,0,0,0,0,0,0,0 

    # Center of mass kinematics     
    if(n_sol==2):
        # Ejectile lab frame 
        P3=momentum_rel(E3,m3)                                      # Ejectile momentum in lab frame
        #print("P3",P3,"E3",E3)

        # Recoil lab frame 
        E4=Etot-E3-m3-m4_new                                        # Recoil energy in lab frame
        P4=math.sqrt(math.pow(E4+m4_new,2)-math.pow(m4_new,2))      # Recoil momentum in lab frame
        #print("m1",m1,"E1",E1,"m3",m3,"E3",E3)

        # somtimes there's a rounding error when cos_arg is exactly 1 (ie theta4=0)
        cos_arg = (m1*E1 + m4_new*E4 - m3*E3) / (2 * math.sqrt(m1 * m4_new * E1 * E4))
        if(cos_arg<1):  # Select "good events": ie physical values for theta4
            theta4=math.acos(cos_arg)   # Recoil theta angle in lab frame (radian)
        elif(1<cos_arg<1.000001):   # rounding error
            theta4=0
        else:   # if truly unphysical, throw error
            theta4=0
            thetaOK=0
        
        #print("E4",E4,"P4",P4,"theta4",theta4*180/math.pi)
        # '''''''''' Center of mass frame  ''''''''''''
        # Velocities
        Beta_CM_ini=P1/Etot                                         # Center of mass initial velocity - beta 
        E_CM=E1*m2/(m1+m2)                                          # Center of mass energy 
        Gamma_CM_ini=1./math.sqrt(1-math.pow(Beta_CM_ini,2))        # Center of mass initial velocity - gamma         
        Beta_CM_fin=math.sqrt(math.pow((Gamma_CM_ini*(m1+m2)*Beta_CM_ini)/(m3+m4_new),2)/(1+math.pow((Gamma_CM_ini*(m1+m2)*Beta_CM_ini)/(m3+m4_new),2)))     # Center of mass initial velocity - beta
        Gamma_CM_fin=1./math.sqrt(1-math.pow(Beta_CM_fin,2))        # Center of mass final velocity - gamma            
        #print("Beta_CM_ini",Beta_CM_ini,"Gamma_CM_ini",Gamma_CM_ini)         
        #print("Beta_CM_fin",Beta_CM_fin,"Gamma_CM_fin",Gamma_CM_fin)
        # Energies
        E3_CM=Gamma_CM_fin*(E3+m3-P3*math.cos(theta3)*Beta_CM_fin)      # Ejectile CM energy
        P3_CM=math.sqrt(math.pow(E3_CM,2)-math.pow(m3,2))               # Ejectile CM momentum
        E4_CM=Gamma_CM_fin*(E4+m4_new-P4*math.cos(theta4)*Beta_CM_fin)  # Recoil CM energy
        P4_CM=math.sqrt(math.pow(E4_CM,2)-math.pow(m4_new,2))           # Recoil CM momentum
        #print("E3CM",E3_CM,"P3_CM",P3_CM,"E4_CM",E4CM,"P4_CM",P4_CM)
        # Angles
        theta3_CM=math.atan((P3*math.sin(theta3))/(Gamma_CM_fin*(P3*math.cos(theta3)-(E3+m3)*Beta_CM_fin)))*180/math.pi # Ejectile CM theta angle
        if(theta3_CM<0):
            theta3_CM=theta3_CM+180
        theta4_CM=180-theta3_CM                                         # Recoil CM theta angle        
        #print("theta3_CM",theta3_CM,"theta4_CM",theta4_CM)
        phi3_lab = np.random.uniform(0.0, math.pi)                      # Ejectile lab phi angle
        phi4_lab = -phi3_lab;                                           # Recoil lab phi angle
        if(phi4_lab<0):
            phi4_lab = phi4_lab+2*math.pi
        # Momenta
        pX_3,pY_3,pZ_3=SetMagThetaPhi(P3,theta3,phi3_lab)              # Ejectile lab momentum x,y,z
        pX_4,pY_4,pZ_4=SetMagThetaPhi(P4,theta4,phi4_lab)              # Recoil lab momentum x,y,z
        #print("pX_4",pX_4,"pY_4",pY_4,"pZ_4",pZ_4)
        #print("pX_3",pX_3,"pY_3",pY_3,"pZ_3",pZ_3)
        # Checking variables  
        #P3_check=Mag(pX_3,pY_3,pZ_3)
        #E3_check=math.sqrt(math.pow(P3_check,2)+math.pow(m3,2))
        #print("P3 after CM",P3_check,"E3 after CM",E3_check-m3)         
        #phi_3,theta_3 = get_theta_phi(pX_3,pY_3,pZ_3)        
        #print("theta3_check",theta3*180/math.pi)    
        #
        #P4_check=Mag(pX_4,pY_4,pZ_4)
        #E4_check=math.sqrt(math.pow(P4_check,2)+math.pow(m4_new,2))
        #print("P4 after CM",P4_check,"E4 after CM",E4_check-m4_new)         
        #phi_4,theta_4 = get_theta_phi(pX_4,pY_4,pZ_4)
        #print("theta4_check",theta_4*180/math.pi)

        if(thetaOK==0):
            print("test failed: recoil angle unphysical, cos argument: {:.8f}".format(cos_arg))
            E3,E4,pX_4,pY_4,pZ_4,pX_3,pY_3,pZ_3,theta4,theta3,theta4_CM,theta3_CM,phi4_lab=0,0,0,0,0,0,0,0,0,0,0,0,0

    return E3,E4,pX_4,pY_4,pZ_4,pX_3,pY_3,pZ_3,theta4,theta3,theta4_CM,theta3_CM,phi4_lab,thetaOK

# _______________________________________________________________________________________________
# GAMMA DEEXCITATION FUNCTION : A relativistic version will have to be done at some point

def dexec_gamma(mass_recoil,q_recoil,Ek_recoil,Eex_before,E_gamma):

    if E_gamma<=0 : 
        E_gamma=Eex_before #desexcitation gamma

    # CM velocity
    V_CM = math.sqrt(2*mass_recoil*Ek_recoil)/mass_recoil   
    #print("V_CM",V_CM)
    # Heavy residue energy (after gamma emission) in the CM frame
    E_HR_CM = math.pow(E_gamma,2.)/(2.*mass_recoil) 
    #print("E_HR_CM",E_HR_CM)           
    # Heavy residue momentum (after gamma emission) in the CM frame
    P_HR_CM=math.sqrt(2*E_HR_CM*mass_recoil)
    #print("P_HR_CM",P_HR_CM)
   
    # Excitation energy after gamma emission
    Eex_after = Eex_before-E_HR_CM-E_gamma 

    # Sampling on theta and phi in CM
    theta_HR_CM = math.acos(np.random.uniform(-1.,1.))                      # CM theta 
    phi_HR_CM = np.random.uniform(0.,2*math.pi)                             # CM phi

    pX_HR,pY_HR,pZ_HR=SetMagThetaPhi(P_HR_CM,theta_HR_CM,phi_HR_CM)         # HR CM momentum 3 axis         

    # Back in the LAB frame    
    pZ_HR = pZ_HR+mass_recoil*V_CM                                          # Adding CM velocity to pz component      
    P_HR=Mag(pX_HR,pY_HR,pZ_HR) 
    #print("P_HR",P_HR)                                                      # LAB momentum
    Ek_HR=math.pow(P_HR,2.)/(2.*mass_recoil)                                # LAB kinetic energy
    #print("Ek_HR",Ek_HR)
    phi_LAB,theta_LAB=get_theta_phi(pX_HR,pY_HR,pZ_HR)                      # LAB theta and phi
    #print("phi_LAB",phi_LAB)
    #print("theta_LAB",theta_LAB*180/math.pi)
    #print("Ek_HR ",Ek_HR)
    
    return Ek_HR

# _______________________________________________________________________________________________
# NEUTRON DEEXCITATION FUNCTION : A relativistic version will have to be done at some point

def dexec_neutron(mass_recoil,q_recoil,A_recoil,Ek_recoil,Eex,Sn,mass_neutron):

    Eex_before = Eex    
    En = -1  

    # Heavy residue 
    A_HR = A_recoil-1
    Z_HR = q_recoil
    mass_HR = atomic_to_nuclear_mass(nuclide_data.weight(Z_HR,A_HR), Z_HR)*amu

    # CM velocity
    V_CM = math.sqrt(2*mass_recoil*Ek_recoil)/mass_recoil 

    # Maximum energy of the neutron
    En_max = mass_HR/(mass_HR+mass_neutron)*(Eex_before-Sn) #maximum energy of the neutron
        
    # to be improved
    # not an uniform distribution!! Maxwel with a given temperature depending on E*
    if(En<=0):
        En = np.random.uniform(0.025,En_max)

    # Excitation energy after gamma emission
    Eex_after = Eex_before-Sn-En*(1.+mass_neutron/mass_HR)
    E_HR_CM = mass_neutron/mass_HR*En
    P_HR_CM = math.sqrt(2.*mass_HR*E_HR_CM)

    # Sampling on theta and phi in CM
    theta_HR_CM = math.acos(np.random.uniform(-1.,1.))                      # CM theta 
    phi_HR_CM = np.random.uniform(0,2*math.pi)                              # CM phi
    
    pX_HR,pY_HR,pZ_HR=SetMagThetaPhi(P_HR_CM,theta_HR_CM,phi_HR_CM)         # HR CM momentum 3 axis 
    pZ_HR2=mass_recoil*V_CM 
    # Back in the LAB frame    
    pZ_HR = pZ_HR+mass_HR*V_CM                                              # Adding CM velocity to pz component
       
    P_HR=Mag(pX_HR,pY_HR,pZ_HR)                                             # LAB momentum
    Ek_HR=math.pow(P_HR,2)/(2*mass_HR)                                      # LAB kinetic energy
    
    phi_LAB,theta_LAB=get_theta_phi(pX_HR,pY_HR,pZ_HR)
    #print("theta_LAB n",theta_LAB*180/math.pi)

    return Ek_HR,mass_HR,Eex_after,pZ_HR