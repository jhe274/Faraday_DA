import numpy as np
from dataclasses import dataclass, field

@dataclass(frozen=True)
class Constants:
    """
    A class containing universal constants, Potassium isotopes' properties,
    and electronic g-factors.
    """
    # Universal Constants
    r_e: float = 2.8179403227e-15         # Classical electron radius [m]
    mu_B: float = 9.2740100783e-24        # Bohr magneton [J/T]
    h: float = 6.62607015e-34             # Planck constant [J/Hz]
    k_B: float = 1.380649e-23             # Boltzmann constant [J/K]
    c: int = 299792458                    # Speed of light [m/s]
    e: float = 1.602176634e-19            # Elementary charge [C]
    m_K39: float = 1.6605402e-27 * 38.96370668  # Potassium-39 atomic mass [kg]
    m_e: float = 9.1093837015e-31         # Electron mass [kg]
    alpha: float = r_e * c / 6            # Global prefactor

    # Potassium-39 properties (D1/D2 lines)
    K39_D1_Hz: float = 389286.058716e9   # [Hz] D1 line frequency in vacuum
    K39_D2_Hz: float = 391016.17003e9    # [Hz] D2 line frequency in vacuum
    K39_D1_m: float = c / K39_D1_Hz      # [m] D1 line wavelength in vacuum
    K39_D2_m: float = c / K39_D2_Hz      # [m] D2 line wavelength in vacuum
    K39_D2_A_Hz: float = 391015.99413e9  # [Hz] |F=2> -> |F'=3,2,1> in vacuum
    K39_D2_C_Hz: float = 391016.44456e9  # [Hz] |F=1> -> |F'=2,1,0> in vacuum
    K39_D2_B_Hz: float = (K39_D2_A_Hz + K39_D2_C_Hz) / 2  # Ground state crossover [Hz]

    # Potassium-41 properties (D1/D2 lines)
    K41_D1_Hz: float = 389286.294205e9   # [Hz] D1 line frequency in vacuum
    K41_D2_Hz: float = 391016.40621e9    # [Hz] D2 line frequency in vacuum
    K41_D1_m: float = c / K41_D1_Hz      # [m] D1 line wavelength in vacuum
    K41_D2_m: float = c / K41_D2_Hz      # [m] D2 line wavelength in vacuum
    K41_D2_A_Hz: float = 391016.54544e9  # [Hz] |F=2> -> |F'=3,2,1> in vacuum
    K41_D2_C_Hz: float = 391016.79394e9  # [Hz] |F=1> -> |F'=2,1,0> in vacuum
    K41_D2_B_Hz: float = (K41_D2_A_Hz + K41_D2_C_Hz) / 2  # Ground state crossover [Hz]

    # Electronic g-factors
    g_factors: dict = field(default_factory=lambda: {
        "D1": 2/3, 
        "D2": 4/3, 
        "ground": 2.00229421})