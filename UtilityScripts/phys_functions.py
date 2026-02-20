#!/usr/bin/env python3.5
import math
import numpy as np
import scipy.constants as const
import os


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

#GLOBAL variables
m0 = const.value(u'atomic mass constant energy equivalent in MeV');     # unit mass in MeV/c^2
e0 = const.e
c = const.c

def Q_val_proj_targ(m_proj, m_targ, m_final):       # Q value for projectile + target -> final nucleus
    return m_proj + m_targ - m_final

def Q_val_proj_targ_ejec(m_proj, m_targ, m_eject, m_final):       # Q value for projectile + target -> ejectile + final nucleus
    return  m_proj + m_targ - m_eject - m_final

def momentum(Ek, m):    # non-relativistic momentum in MeV c
    return math.sqrt(2 * m * Ek)

def momentum_rel(Ek, m):
    if Ek <= 0.0:       # if kinetic energy is negative, return 0
        print("Kinetic energy is negative: Ek = ", Ek)
        return 0.0
    return math.sqrt(Ek * (Ek + 2.0 * m))

def GetBetaVec(FV):
    """
    Return beta vector (v/c) from a four-vector FV = (E, px, py, pz)
    """
    E = FV[0]
    p = FV[1:4]

    if E == 0:
        raise ValueError("Energy is zero, cannot compute beta.")

    return p / E


def GetGamma(FV):
    """
    Return Lorentz gamma from a four-vector FV = (E, px, py, pz)
    """
    E = FV[0]
    p = FV[1:4]

    p2 = np.dot(p, p)
    m2 = E*E - p2

    if m2 <= 0:
        raise ValueError("Invariant mass squared <= 0 (check units or FV).")

    m = math.sqrt(m2)
    return E / m

def lorentz_boost(FV, beta_vec, gamma):     # lorentz boost from CM to LAB frame
    if np.linalg.norm(beta_vec) == 0:       # if beta vector is zero, return the four-vector unchanged
        return FV
    
    E_CM, p_CM = FV[0], FV[1:4]             # extract energy and momentum from four-vector
    beta = np.linalg.norm(beta_vec)         # beta vector magnitude
    beta_hat = beta_vec / beta              # beta vector unit vector
    
    p_par = np.dot(p_CM, beta_hat) * beta_hat    # parallel component of the momentum
    p_perp = p_CM - p_par                        # perpendicular component of the momentum
    
    # Lorentz transformation
    E_lab = gamma * (E_CM + np.dot(p_CM, beta_vec))    # energy in the LAB frame
    p_par_lab = gamma * (p_par + beta_vec * E_CM)    # parallel component of the momentum in the LAB frame
    p_lab = p_perp + p_par_lab              # total momentum in the LAB frame
    
    return np.concatenate(([E_lab], p_lab))

def energy_rel(p, A, m0):    # relativistic energy in MeV
    return math.sqrt(p**2 + (A * m0)**2) - (A * m0)

def brho(p, q):    # return brho in T m
    return p * e0 / c * 10**6 / (q * e0)    # return brho in T m

def x_alpha(max_X, max_A):
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
    
def spherical_to_cartesian(mag, theta, phi):
    amag = abs(mag)
    fX = amag * math.sin(theta) * math.cos(phi)
    fY = amag * math.sin(theta) * math.sin(phi)
    fZ = amag * math.cos(theta)
    return np.array([fX, fY, fZ])

def magnitude(fx, fy, fz):
    return math.sqrt(fx**2 + fy**2 + fz**2)


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

# Electron binding energies cache (loaded once)
_binding_energies = None

def _load_binding_energies():       # load electron binding energies from data file (cached)
    global _binding_energies
    if _binding_energies is not None:
        return _binding_energies
    
    _binding_energies = {}
    # File is now in NuclideDataMaster directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    binding_energy_file = os.path.join(script_dir, "..", "NuclideDataMaster", "Electron_Binding_Energies-2007.dat")
    
    if not os.path.exists(binding_energy_file):
        # Return empty dict if file not found - function will handle gracefully
        return _binding_energies
    
    with open(binding_energy_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 2:
                    Z = int(parts[0])
                    total_binding_energy = float(parts[1])  # in eV
                    _binding_energies[Z] = total_binding_energy
    
    return _binding_energies

def atomic_to_nuclear_mass(atomic_mass_amu, Z):     # Convert atomic mass to nuclear mass using electron binding energies.
    binding_energies = _load_binding_energies()
    
    if Z not in binding_energies:
        # For Z=1 and Z=2, return atomic mass (no correction available)
        # This is a reasonable approximation since electron binding energy is very small for light elements
        return atomic_mass_amu
    
    # Get electron binding energy in eV
    binding_energy_ev = binding_energies[Z]
    
    # Convert eV to amu (1 eV = 1.07354410233e-9 amu)
    # Using E = mc^2 with c = 299792458 m/s and 1 amu = 931.49410242 MeV/c^2
    binding_energy_amu = binding_energy_ev * 1.07354410233e-9
    
    # Nuclear mass = atomic mass - electron binding energy
    nuclear_mass_amu = atomic_mass_amu - binding_energy_amu
    
    return nuclear_mass_amu
