#!/usr/bin/env python
# --
# Nuclear Mass Converter
# Converts atomic masses to nuclear masses using electron binding energies
# Last modification : G. Leckenby - 8/07/2025

from __future__ import print_function
import numpy as np
import os

class NuclearMassConverter(object):
    """
    Class to convert atomic masses to nuclear masses using electron binding energies.
    Uses the Electron_Binding_Energeies-2007.dat file for binding energy data.
    """
    
    def __init__(self, binding_energy_file="Electron_Binding_Energeies-2007.dat"):
        """
        Initialize the converter with electron binding energy data.
        
        Parameters:
        -----------
        binding_energy_file : str
            Path to the electron binding energies data file
        """
        self.binding_energy_file = binding_energy_file
        self.binding_energies = {}
        self.load_binding_energies()
    
    def load_binding_energies(self):
        """
        Load electron binding energies from the data file.
        """
        if not os.path.exists(self.binding_energy_file):
            raise IOError("Binding energy file not found: {}".format(self.binding_energy_file))
        
        with open(self.binding_energy_file, 'r') as f:
            lines = f.readlines()
        
        # Skip comment lines and header
        for line in lines:
            if line.strip() and not line.startswith('#'):
                # Parse the line: Z, TotBEn, 1e-ion, 2e-ion, ...
                parts = line.strip().split()
                if len(parts) >= 2:
                    Z = int(parts[0])
                    total_binding_energy = float(parts[1])  # in eV
                    self.binding_energies[Z] = total_binding_energy
    
    def atomic_to_nuclear_mass(self, atomic_mass_amu, Z):
        """
        Convert atomic mass to nuclear mass.
        
        Parameters:
        -----------
        atomic_mass_amu : float
            Atomic mass in atomic mass units (amu)
        Z : int
            Atomic number (number of protons)
        
        Returns:
        --------
        float
            Nuclear mass in atomic mass units (amu)
        """
        if Z not in self.binding_energies:
            # For Z=1 and Z=2, return atomic mass (no correction available)
            # This is a reasonable approximation since electron binding energy is very small for light elements
            return atomic_mass_amu
        
        # Get electron binding energy in eV
        binding_energy_ev = self.binding_energies[Z]
        
        # Convert eV to amu (1 eV = 1.07354410233e-9 amu)
        # Using E = mc^2 with c = 299792458 m/s and 1 amu = 931.49410242 MeV/c^2
        binding_energy_amu = binding_energy_ev * 1.07354410233e-9
        
        # Nuclear mass = atomic mass - electron binding energy
        nuclear_mass_amu = atomic_mass_amu - binding_energy_amu
        
        return nuclear_mass_amu
    
    def nuclear_to_atomic_mass(self, nuclear_mass_amu, Z):
        """
        Convert nuclear mass to atomic mass.
        
        Parameters:
        -----------
        nuclear_mass_amu : float
            Nuclear mass in atomic mass units (amu)
        Z : int
            Atomic number (number of protons)
        
        Returns:
        --------
        float
            Atomic mass in atomic mass units (amu)
        """
        if Z not in self.binding_energies:
            # For Z=1 and Z=2, return nuclear mass (no correction available)
            # This is a reasonable approximation since electron binding energy is very small for light elements
            return nuclear_mass_amu
        
        # Get electron binding energy in eV
        binding_energy_ev = self.binding_energies[Z]
        
        # Convert eV to amu
        binding_energy_amu = binding_energy_ev * 1.07354410233e-9
        
        # Atomic mass = nuclear mass + electron binding energy
        atomic_mass_amu = nuclear_mass_amu + binding_energy_amu
        
        return atomic_mass_amu

# Global converter instance
_nuclear_mass_converter = None

def get_nuclear_mass_converter():
    """
    Get or create a global nuclear mass converter instance.
    
    Returns:
    --------
    NuclearMassConverter
        Global converter instance
    """
    global _nuclear_mass_converter
    if _nuclear_mass_converter is None:
        _nuclear_mass_converter = NuclearMassConverter()
    return _nuclear_mass_converter

def atomic_to_nuclear_mass(atomic_mass_amu, Z):
    """
    Convenience function to convert atomic mass to nuclear mass.
    
    Parameters:
    -----------
    atomic_mass_amu : float
        Atomic mass in atomic mass units (amu)
    Z : int
        Atomic number (number of protons)
    
    Returns:
    --------
    float
        Nuclear mass in atomic mass units (amu)
    """
    converter = get_nuclear_mass_converter()
    return converter.atomic_to_nuclear_mass(atomic_mass_amu, Z)

def nuclear_to_atomic_mass(nuclear_mass_amu, Z):
    """
    Convenience function to convert nuclear mass to atomic mass.
    
    Parameters:
    -----------
    nuclear_mass_amu : float
        Nuclear mass in atomic mass units (amu)
    Z : int
        Atomic number (number of protons)
    
    Returns:
    --------
    float
        Atomic mass in atomic mass units (amu)
    """
    converter = get_nuclear_mass_converter()
    return converter.nuclear_to_atomic_mass(nuclear_mass_amu, Z)

# Example usage and testing
if __name__ == "__main__":
    # Test the converter
    converter = NuclearMassConverter()
    
    # Example: Convert atomic mass of carbon (Z=6) to nuclear mass
    carbon_atomic_mass = 12.0  # amu
    carbon_Z = 6
    
    carbon_nuclear_mass = converter.atomic_to_nuclear_mass(carbon_atomic_mass, carbon_Z)
    print("Carbon (Z={}):".format(carbon_Z))
    print("  Atomic mass: {} amu".format(carbon_atomic_mass))
    print("  Nuclear mass: {:.12f} amu".format(carbon_nuclear_mass))
    print("  Difference: {:.6f} micro-amu".format((carbon_atomic_mass - carbon_nuclear_mass)*1e6))
    
    # Test reverse conversion
    carbon_atomic_mass_back = converter.nuclear_to_atomic_mass(carbon_nuclear_mass, carbon_Z)
    print("  Atomic mass (recovered): {:.12f} amu".format(carbon_atomic_mass_back))
    print("  Recovery error: {:.2e} amu".format(abs(carbon_atomic_mass - carbon_atomic_mass_back))) 