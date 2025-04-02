
import sys, os
import numpy as np

# This theory calculator is used to plot the BAO points and also to plot the difference in the SN points introduced by the DDR in Remake
# The values in the run should be adjusted accordingly between the two files

# #Calling CAMB
# sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
# camb_installation_path = '/Users/SM1WG/Documents/Research/Cosmology/Codes/Theory/CAMB/'
# camb_path = os.path.realpath(os.path.join(os.getcwd(),camb_installation_path))
# sys.path.insert(0,camb_path)

#For CLASS you just need to have the right class folder compiled in your environment

from classy import Class

c=299792.458

class TheoryCalculator:
    def __init__(self, H0=67.97, ombh2=0.02277699188, omch2=0.1190083805, tau=0.0544, ns=0.9649, As=3.044, omk=0,a0=0.049,a1=-0.0,eta=0,zstar=0.9):
        # Setting up parameters for CLASS
        self.params = {
            'output': 'mPk',  # Optional, depending on which outputs you need (mPk, dTk, etc.)
            'H0': H0,
            'omega_b': ombh2,
            'omega_cdm': omch2,
            'tau_reio': tau,
            'n_s': ns,
            'A_s': np.exp(As) * (10**-10),  # Rescale As if needed
            'alpha0_ddr': a0,
            'alpha1_ddr': a1,
            'zstar_ddr': zstar,
            'eta_function':eta,
            'Omega_k': omk
        }
        
        # Initialize the CLASS instance
        self.cosmo = Class()
        self.cosmo.set(self.params)
        self.cosmo.compute()
    
    def get_all(self):
        # Returns the CLASS cosmology instance
        return self.cosmo

    def get_luminosity_distance(self, z):
        # Luminosity distance in Mpc
        return self.cosmo.angular_distance(z)*(1.+z)*(1.+z)
    
    def get_distance_moduli(self, z):
        # Calculate the distance modulus
        return 5 * np.log10(self.get_luminosity_distance(z)) + 25
    
    def get_luminosity_distance_ddr(self, z):
        # Luminosity distance in Mpc
        return self.cosmo.angular_distance_ddr(z)*(1.+z)*(1.+z)
    
    def get_distance_moduli_ddr(self, z):
        # Calculate the distance modulus
        return 5 * np.log10(self.get_luminosity_distance_ddr(z)) + 25
    
    def get_H(self, z):
        # Get the Hubble parameter at redshift z
        return self.cosmo.Hubble(z) * c  # Convert H(z) from 1/Mpc to km/s/Mpc

    def get_r_drag(self):
        # Sound horizon at drag epoch
        return self.cosmo.rs_drag()
    
    def __del__(self):
        # Clean up the CLASS instance
        self.cosmo.struct_cleanup()
        self.cosmo.empty()


class TheoryCalculatorwp:
    def __init__(self, H0=67.4, ombh2=0.02239657558, omch2=0.1220631924, tau=0.0544, ns=0.9649, As=3.044, omk=0,a0=-0.0,a1=-0.0,eta=0,zstar=0.9,w0=-.851,wa=-0.7):
        # Setting up parameters for CLASS
        self.params = {
            'output': 'mPk',  # Optional, depending on which outputs you need (mPk, dTk, etc.)
            'H0': H0,
            'omega_b': ombh2,
            'omega_cdm': omch2,
            'Omega_Lambda':0,
            'w0_fld':w0,
            'wa_fld':wa,
            'tau_reio': tau,
            'n_s': ns,
            'A_s': np.exp(As) * (10**-10),  # Rescale As if needed
            'alpha0_ddr': a0,
            'alpha1_ddr': a1,
            'zstar_ddr': zstar,
            'eta_function':eta,
            'Omega_k': omk
        }
        
        # Initialize the CLASS instance
        self.cosmo = Class()
        self.cosmo.set(self.params)
        self.cosmo.compute()
    
    def get_all(self):
        # Returns the CLASS cosmology instance
        return self.cosmo

    def get_luminosity_distance(self, z):
        # Luminosity distance in Mpc
        return self.cosmo.angular_distance(z)*(1.+z)*(1.+z)
    
    def get_distance_moduli(self, z):
        # Calculate the distance modulus
        return 5 * np.log10(self.get_luminosity_distance(z)) + 25
    
    def get_luminosity_distance_ddr(self, z):
        # Luminosity distance in Mpc
        return self.cosmo.angular_distance_ddr(z)*(1.+z)*(1.+z)
    
    def get_distance_moduli_ddr(self, z):
        # Calculate the distance modulus
        return 5 * np.log10(self.get_luminosity_distance_ddr(z)) + 25
    
    def get_H(self, z):
        # Get the Hubble parameter at redshift z
        return self.cosmo.Hubble(z) * c  # Convert H(z) from 1/Mpc to km/s/Mpc

    def get_r_drag(self):
        # Sound horizon at drag epoch
        return self.cosmo.rs_drag()
    
    def __del__(self):
        # Clean up the CLASS instance
        self.cosmo.struct_cleanup()
        self.cosmo.empty()

class TheoryCalculatorw:
    def __init__(self, H0=68.03, ombh2=0.02281722206, omch2=0.1199590737, tau=0.0544, ns=0.9649, As=3.044, omk=0,a0=-0.0,a1=-0.0,eta=0,zstar=0.9,w0=-0.827,wa=-0.75):
        # Setting up parameters for CLASS
        self.params = {
            'output': 'mPk',  # Optional, depending on which outputs you need (mPk, dTk, etc.)
            'H0': H0,
            'omega_b': ombh2,
            'omega_cdm': omch2,
            'Omega_Lambda':0,
            'w0_fld':w0,
            'wa_fld':wa,
            'tau_reio': tau,
            'n_s': ns,
            'A_s': np.exp(As) * (10**-10),  # Rescale As if needed
            'alpha0_ddr': a0,
            'alpha1_ddr': a1,
            'zstar_ddr': zstar,
            'eta_function':eta,
            'Omega_k': omk
        }
        
        # Initialize the CLASS instance
        self.cosmo = Class()
        self.cosmo.set(self.params)
        self.cosmo.compute()
    
    def get_all(self):
        # Returns the CLASS cosmology instance
        return self.cosmo

    def get_luminosity_distance(self, z):
        # Luminosity distance in Mpc
        return self.cosmo.angular_distance(z)*(1.+z)*(1.+z)
    
    def get_distance_moduli(self, z):
        # Calculate the distance modulus
        return 5 * np.log10(self.get_luminosity_distance(z)) + 25
    
    def get_luminosity_distance_ddr(self, z):
        # Luminosity distance in Mpc
        return self.cosmo.angular_distance_ddr(z)*(1.+z)*(1.+z)
    
    def get_distance_moduli_ddr(self, z):
        # Calculate the distance modulus
        return 5 * np.log10(self.get_luminosity_distance_ddr(z)) + 25
    
    def get_H(self, z):
        # Get the Hubble parameter at redshift z
        return self.cosmo.Hubble(z) * c  # Convert H(z) from 1/Mpc to km/s/Mpc

    def get_r_drag(self):
        # Sound horizon at drag epoch
        return self.cosmo.rs_drag()
    
    def __del__(self):
        # Clean up the CLASS instance
        self.cosmo.struct_cleanup()
        self.cosmo.empty()