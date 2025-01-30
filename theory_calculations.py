import numpy as np
import math
from constants import Constants as Consts

class Theory:
    
    def __init__(self):
        """
        Initializes the Theory class with physical constants.
        """
        self.consts = Consts()

    def Kn_density(self, T):
        """
        Calculates Potassium number density based on an empirical formula.
        
        Parameters:
            T (float): Temperature in Celsius.
        
        Returns:
            float: Potassium number density [m^-3].
        """
        if 24.85 < T <= 63.35:
            Kp = 10 ** (9.967 - 4646 / (273.15 + T))  # [Pa]
        elif T > 63.35:
            Kp = 10 ** (9.408 - 4453 / (273.15 + T))  # [Pa]
        else:
            raise ValueError("Temperature T must be above 24.85°C")
        
        return Kp / (self.consts.k_B * (273.15 + T))  # [m^-3]
    
    def doppler_broad(self, nu0, T):
        """
        Calculates HWHM Doppler broadening for a given transition.
        
        Parameters:
            nu0 (float): Central frequency [Hz].
            T (float): Temperature in Celsius.
        
        Returns:
            float: Doppler broadening [Hz].
        """
        return nu0 * np.sqrt(
            (8 * self.consts.k_B * (273.15 + T) * math.log(2)) / (self.consts.m_K39 * self.consts.c**2)
        )

    def Zeeman_splitting(self, B, state):
        """
        General Zeeman splitting calculation for any state.

        Parameters:
            B (float): Magnetic field strength [T].
            state (str): The atomic state ('D1', 'D2', or 'ground').

        Returns:
            float: Zeeman energy shift [J].
        """
        if state not in self.consts.g_factors:
            raise ValueError("Invalid state. Choose from 'D1', 'D2', or 'ground'.")
        
        g_factor = self.consts.g_factors[state]
        return g_factor * 0.5 * self.consts.mu_B * B  # [J]
    
    def diamagnetic_FR(self, nu, l, Kn, T, B, gamma_D1, gamma_D2):
        """
        Calculates the Diamagnetic Faraday Rotation.
        
        Parameters:
            nu (float): Probe frequency [Hz].
            l (float): Path length [m].
            Kn (float): Number density [m^-3].
            T (float): Temperature [°C].
            B (float): Magnetic field strength [T].
            gamma_D1, gamma_D2 (float): Linewidths [MHz].
        
        Returns:
            float: Rotation angle [rad].
        """
        delta_nu_D1 = nu - self.consts.K39_D1_Hz
        delta_nu_D2 = nu - self.consts.K39_D2_Hz

        dia_FR = ((7 * (delta_nu_D2**2 - (gamma_D2 * 1e6)**2 / 4) / (delta_nu_D2**2 + (gamma_D2 * 1e6)**2 / 4)**2) + 
                  (4 * (delta_nu_D1**2 - (gamma_D1 * 1e6)**2 / 4) / (delta_nu_D1**2 + (gamma_D1 * 1e6)**2 / 4)**2)) / (3 * self.consts.h)

        return self.consts.alpha * Kn * 1e14 * l * self.consts.mu_B * np.sign(B) * B * 1e-4 * dia_FR  # [rad]
    
    def paramagnetic_FR(self, nu, l, Kn, T, P, gamma_D1, gamma_D2):
        """
        Calculates the Paramagnetic Faraday Rotation.
        
        Parameters:
            nu (float): Probe frequency [Hz].
            l (float): Path length [m].
            Kn (float): Number density [m^-3].
            T (float): Temperature [°C].
            P (float): Polarization parameter.
            gamma_D1, gamma_D2 (float): Linewidths [MHz].
        
        Returns:
            float: Rotation angle [rad].
        """
        delta_nu_D1 = nu - self.consts.K39_D1_Hz
        delta_nu_D2 = nu - self.consts.K39_D2_Hz

        para_FR = P * ((delta_nu_D2 / ((delta_nu_D2 - self.doppler_broad(self.consts.K39_D2_Hz, T))**2 + (gamma_D2 * 1e6)**2 / 4)) - 
                       (delta_nu_D1 / ((delta_nu_D1 - self.doppler_broad(self.consts.K39_D1_Hz, T))**2 + (gamma_D1 * 1e6)**2 / 4)))
        
        return self.consts.alpha * Kn * 1e14 * l * para_FR  # [rad]
    
    def resonant_FR(self, nu, Kn, T, B, P, gamma_D1, gamma_D2, const):
        """
        Calculates the Resonant Faraday Rotation.
        
        Parameters:
            nu (float): Probe frequency [Hz].
            Kn (float): Number density [m^-3].
            T (float): Temperature [°C].
            B (float): Magnetic field strength [T].
            P (float): Polarization parameter.
            gamma_D1, gamma_D2 (float): Linewidths [MHz].
            const (float): Additional correction factor.
        
        Returns:
            float: Rotation angle [rad].
        """
        l = 7.182e-2  # [m] Path length of the vapor cell used in the experiment
        return (self.diamagnetic_FR(nu, l, Kn, T, B, gamma_D1, gamma_D2) + 
                self.paramagnetic_FR(nu, l, Kn, T, P, gamma_D1, gamma_D2) + const)