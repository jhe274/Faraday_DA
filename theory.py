import numpy as np
import math
from constants import Constants as Consts

class Theory:
    
    def __init__(self):
        self.consts = Consts()

    def Kn_density(self, T):
        '''
        Potassium number density based on experienced formula
        '''
        if 24.85 < T <= 63.35:
            Kp = 10 ** (9.967 - 4646 / (273.15 + T))                                                                                    # [Pa]
        elif 63.35 < T:
            Kp = 10 ** (9.408 - 4453 / (273.15 + T))                                                                                    # [Pa]
        else:
            raise ValueError("Temperature T must be above 24.85°C")
        
        return Kp / (self.consts.k_B * (273.15 + T))                                                                                    # [m^-3]
    
    def doppler_broad(self, nu0, T):
        Delta_D = nu0 * np.sqrt((8 * self.consts.k_B * (273.15+T) * math.log(2)) / (self.consts.m_K39 * self.consts.c**2))

        return Delta_D

    def Zeeman_ground(self, B):
        '''
        Calculations of Zeeman splitting 
        '''
        return self.consts.g_G * (1/2) * self.consts.mu_B * B

    def Zeeman_D1(self, B):
        """
        Zeeman splitting at potassium D1 line
        """
        E_D1 = self.consts.g_D1 * (1/2) * self.consts.mu_B * B
        Delta_E_D1 = -self.Zeeman_ground(B) - E_D1 - (self.Zeeman_ground(B) + E_D1)
        Delta_Lambda_D1 = - pow(self.consts.Lambda39_D1,2) * Delta_E_D1 / (self.consts.h * self.consts.c)                               # [m]
        Delta_nu_D1 = self.consts.Nu39_D1 - (self.consts.c / (self.consts.Lambda39_D1 + Delta_Lambda_D1/2))*2                           # [Hz]

        return Delta_Lambda_D1, Delta_nu_D1

    def Zeeman_D2(self, B):
        """
        Zeeman splitting at potassium D2 line
        """
        E_D2 = self.consts.g_D2 * (3/2) * self.consts.mu_B * B
        Delta_E_D2 = -self.Zeeman_ground(B) - E_D2 - (self.Zeeman_ground(B) + E_D2)
        Delta_Lambda_D2 = - pow(self.consts.Lambda39_D2,2) * Delta_E_D2 / (self.consts.h * self.consts.c)                               # [m]
        Delta_nu_D2 = self.consts.Nu39_D2 - (self.consts.c / (self.consts.Lambda39_D2 + Delta_Lambda_D2/2)) * 2                         # [Hz]

        return Delta_Lambda_D2, Delta_nu_D2
    
    def diamagnetic_FR(self, nu, l, Kn, T, B):
        """
        Calculations of diamagnetic Faraday rotation
        """
        Delta_Lambda_D1, Delta_nu_D1 = self.Zeeman_D1(B)
        Delta_Lambda_D2, Delta_nu_D2 = self.Zeeman_D2(B)

        print(f"Zeeman splitting of K D1 line with B={round(B*pow(10,4),3)} G is " , 
              str(Delta_Lambda_D1*pow(10,9)), " nm = ", str(Delta_nu_D1*pow(10,-9)), " GHz")
        print(f"Zeeman splitting of K D2 line with B={round(B*pow(10,4),3)} G is ",
              str(Delta_Lambda_D2*pow(10,9)), " nm = ", str(Delta_nu_D2*pow(10,-9)), " GHz")
        
        delta_nu_D2 = nu - self.consts.Nu39_D2
        delta_nu_D1 = nu - self.consts.Nu39_D1
        delta_doppler_D2 = self.doppler_broad(self.consts.Nu39_D2, T)
        delta_doppler_D1 = self.doppler_broad(self.consts.Nu39_D1, T)

        dia_FR1 = ( 
                (7*(delta_nu_D2**2 - delta_doppler_D2**2/4) / (delta_nu_D2**2 + delta_doppler_D2**2/4)**2) + 
                (4*(delta_nu_D1**2 - delta_doppler_D1**2/4) / (delta_nu_D1**2 + delta_doppler_D1**2/4)**2) - 
                (2 / (delta_nu_D2 * delta_nu_D1))) / (3 * self.consts.h)
        
        # dia_FR2 = np.sign(B) * ((nu / (self.consts.Nu39_D1 * (nu - self.consts.Nu39_D1))) 
            #                            - (nu / (self.consts.Nu39_D2 * (nu - self.consts.Nu39_D2)))) / (self.consts.k_B * (273.15 + T))
        
        dia_theta = self.consts.alpha * Kn * 1e14 * l * self.consts.mu_B * np.sign(B) * B * 1e-4 * dia_FR1                      # [rad]

        return dia_theta
    
    def paramagnetic_FR(self, nu, l, Kn, T, P):
        """
        Calculations of paramagnetic Faraday rotation
        """

        delta_nu_D2 = nu - self.consts.Nu39_D2
        delta_nu_D1 = nu - self.consts.Nu39_D1
        delta_doppler_D2 = self.doppler_broad(self.consts.Nu39_D2, T)
        delta_doppler_D1 = self.doppler_broad(self.consts.Nu39_D1, T)

        para_FR = (delta_nu_D2 / ((delta_nu_D2 - delta_doppler_D2) ** 2 + self.doppler_broad(self.consts.Nu39_D2, T)**2/4))
        - (delta_nu_D1 / ((delta_nu_D1 - delta_doppler_D1) ** 2 + self.doppler_broad(self.consts.Nu39_D1, T)**2/4))
        
        para_theta = self.consts.alpha * Kn * 1e14 * l * P * para_FR                                                            # [rad]

        return para_theta
    
    def resonant_FR(self, nu, Kn, T, B, P, const):
        l = (7.5-0.159*2)*1e-2                                                                                                  # [m]
        theta = self.diamagnetic_FR(nu, l, Kn, T, B) + self.paramagnetic_FR(nu, l, Kn, T, P) + const

        return theta