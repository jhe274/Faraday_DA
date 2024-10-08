import numpy as np

class Constants:
    
    def __init__(self):
        '''
        Universal constants
        '''
        self.r_e = np.float64(2.8179403227 * 1e-15)                                              # Classical electron radius: [m]
        self.mu_B = np.float64(9.2740100783 * 1e-24)                                             # Bohr magneton: [J/T]
        self.h = np.float64(6.62607015 * 1e-34)                                                  # Planck constant: [J/Hz]
        self.k_B = np.float64(1.380649 * 1e-23)                                                  # Boltzmann constant: [J/K]
        self.c = 299792458                                                                       # Speed of light: [m/s]
        self.e = np.float64(1.602176634 * 1e-19)                                                 # Elementary charge: [C]
        self.m_K39 = np.float64(1.6605402*1e-27 * 38.96370668)                                   # Potassium 39 atom mass: [kg]
        self.m_e = np.float64(9.1093837015 * 1e-31)                                              # Electron mass: [kg]
        self.alpha = self.r_e * self.c / 6                                                       # Global prefactor

        '''
        Potassium K39 properties
        '''
        self.Nu39_D1 = 389286.058716 * 1e9                                                       # [Hz], K39 D1 line frequency in vacuum
        self.Nu39_D2 = 391016.17003 * 1e9                                                        # [Hz], K39 D2 line frequency in vacuum
        self.Lambda39_D1 = self.c / (self.Nu39_D1)                                               # [m], K39 D1 line wavelength in vacuum
        self.Lambda39_D2 = self.c / (self.Nu39_D2)                                               # [m], K39 D2 line wavelength in vacuum
        self.Nu39_D2_A = 391015.99413 * 1e9                                                      # [Hz], K39 D2 line, |F=2> -> |F'=3,2,1> in vacuum
        self.Nu39_D2_C = 391016.44456 * 1e9                                                      # [Hz], K39 D2 line, |F=1> -> |F'=2,1,0> in vacuum
        self.Nu39_D2_B = (self.Nu39_D2_A + self.Nu39_D2_C) / 2                                   # [Hz], K39 D2 line, grouns state cross over in vacuum

        '''
        Potassium K41 properties
        '''
        self.Nu41_D1 = 389286.294205 * 1e9                                                       # [Hz], K41 D1 line frequency in vacuum
        self.Nu41_D2 = 391016.40621 * 1e9                                                        # [Hz], K41 D2 line frequency in vacuum
        self.Lambda41_D1 = self.c / (self.Nu41_D1)                                               # [m], K41 D1 line wavelength in vacuum
        self.Lambda41_D2 = self.c / (self.Nu41_D2)                                               # [m], K41 D2 line wavelength in vacuum
        self.Nu41_D2_A = 391016.54544 * 1e9                                                      # [Hz], K41 D2 line, |F=2> -> |F'=3,2,1> in vacuum
        self.Nu41_D2_C = 391016.79394 * 1e9                                                      # [Hz], K41 D2 line, |F=1> -> |F'=2,1,0> in vacuum
        self.Nu41_D2_B = (self.Nu41_D2_A + self.Nu41_D2_C) / 2                                   # [Hz], K41 D2 line, grouns state cross over in vacuum


        '''
        Potassium total electronic g-factors
        '''
        self.g_D1 = 2/3
        self.g_D2 = 4/3
        self.g_G = 2.00229421