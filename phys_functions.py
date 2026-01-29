#!/usr/bin/env python
import os.path
import warnings
import re
import gzip
import copy
import string
from functools import total_ordering
import math
import numpy as np
import uncertainties as unc #this package was installed to be able to use the next line
from nuclidedatamaster import nuclide_data #this requires the nuclidedatamaster directory
import scipy.constants as const


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

#GLOBAL variables
m0=const.value(u'atomic mass constant energy equivalent in MeV'); #unit mass in MeV/c^2
e0=const.e
c=const.c

def Q_val_proj_targ(m_proj,m_targ,m_final):
    m_Q_reac = m_proj+m_targ-m_final
    
    return Q_proj_targ

def Q_val_proj_targ_ejec(m_proj,m_targ,m_eject,m_final):
    m_Q_reac = m_proj+m_targ-m_eject-m_final
    
    return m_Q_reac


#def momentum(Ek,A):    # Energy in MeV, momentum in MeV c
#    gamma=Ek/(A*m0)+1
#    beta=math.sqrt(1-1/gamma**2)

def momentum(Ek,m):    # Energy in MeV, momentum in MeV c
    p=math.sqrt(2*m*Ek)

    return p

def momentum_rel(Ek,m):    # Energy in MeV, momentum in MeV c
    gamma=Ek/m+1
    beta=math.sqrt(1-1/math.pow(gamma,2))
    p=gamma*m*beta

    return p


def energy(p,A,m0):    # Energy in MeV, momentum in MeV c
    
    return math.sqrt(p**2+(A*m0)**2)-(A*m0)


def brho(p,q):    # Energy in MeV, momentum in MeV c
    return p*e0/c*10**6/(q*e0)


def x_alpha(max_X,max_A):
    val=10
    while val>1:
        X=np.random.uniform(-max_X,max_X)
        A=np.random.uniform(-max_A,max_A)
        val=pow(A/max_A,2.)+pow(X/max_X,2.)
    return A,X    

def y_beta(max_Y,max_B):
    val=10
    while val>1:
        Y=np.random.uniform(-max_Y,max_Y)
        B=np.random.uniform(-max_B,max_B)
        val=pow(B/max_B,2.)+pow(Y/max_Y,2.)
        
    return B,Y    
    
def SetMagThetaPhi(mag,theta, phi):
    amag = abs(mag)
    fX = amag*math.sin(theta)*math.cos(phi)
    fY = amag*math.sin(theta)*math.sin(phi)
    fZ = amag*math.cos(theta)
    return fX,fY,fZ


def velocity(P,theta,phi,m):
    vx = (P*math.sin(theta)*math.cos(phi))/m
    vy = P*math.sin(theta)*math.sin(phi)/m
    vz = P*math.cos(theta)/m
    return vx,vy,vz

def Mag(fx,fy,fz):
    rho=math.sqrt(fx*fx+fy*fy+fz*fz)
    return rho

def get_theta_phi(fx,fy,fz):
    rho=Mag(fx,fy,fz)
    phi=math.atan2(fy,fx)
    theta=math.acos(fz/rho)
    return phi,theta


def RotateX(fx,fy,fz,theta):
    x=fx
    y=fy*math.cos(theta)-fz*math.sin(theta)
    z=fy*math.sin(theta)+fz*math.cos(theta)
    return x,y,z

def RotateY(fx,fy,fz,theta):
    x=fx*math.cos(theta)+fz*math.sin(theta)
    y=fy
    z=-fx*math.sin(theta)+fz*math.cos(theta)
    return x,y,z

def RotateZ(fx,fy,fz,theta):
    x=fx*math.cos(theta)-fy*math.sin(theta)
    y=fx*math.sin(theta)+fy*math.cos(theta)
    z=fz
    return x,y,z


def pol1(x, m, b):
    return m*x+b

# Function to convert Z et A as g4beamline understands it (Z_in and A_in are integer values, 
# Z_out and A_out are values converted for g4beamline) 
def define_PDGid(Z_in,A_in):
    Z_out=0
    A_out=0
    
    if Z_in<10:
        Z_out="00"+str(Z_in)
    elif 9<Z_in<100:
        Z_out="0"+str(Z_in)
    else:
        print("ERROR : PDGid definition failed")
        Z_out="000"
        
    if A_in<10:
        A_out="00"+str(A_in)
    elif 9<A_in<100:
        A_out="0"+str(A_in)
    elif 99<A_in<296:
        A_out=""+str(A_in)
    else:
        print("ERROR : PDGid definition failed")
        A_out="000"
        
    return "100"+Z_out+A_out+"0"

