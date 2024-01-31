import numpy as np

class Constants:
    
    def __init__(self):
        '''
        Universal constants
        '''
        self.r_e = np.float64(2.8179403227 * 1e-15)                                             # [m]
        self.mu_B = np.float64(9.2740100783 * 1e-24)                                            # [J/T]
        self.h = np.float64(6.62607015 * 1e-34)                                                 # [J/Hz]
        self.k_B = np.float64(1.380649 * 1e-23)                                                 # [J/K]
        self.c = 299792458                                                                      # [m/s]
        self.alpha = self.r_e * self.mu_B / (6 * self.h)

        '''
        Potassium properties
        '''
        self.Lambda_D1 = 770.108385049 * 1e-9                                                   # [m], K39 D1 line wavelength in vacuum
        self.Lambda_D2 = 766.700921822 * 1e-9                                                   # [m], K39 D2 line wavelength in vacuum
        # self.Lambda_D2 = 766.7005 * 1e-9
        self.Nu_D1 = self.c / self.Lambda_D1                                                    # [Hz], K39 D1 line frequency in vacuum
        self.Nu_D2 = self.c / self.Lambda_D2                                                    # [Hz], K39 D2 line frequency in vacuum
        self.CL_D1 = 770.108385 * 1e-9                                                          # [m], measured center of potassium D1 line
        self.CL_D2 = 766.7009 * 1e-9                                                            # [m], measured center of potassium D2 line

        '''
        Potassium total electronic g-factors
        '''
        self.g_D1 = 2/3
        self.g_D2 = 4/3
        self.g_G = 2.00229421