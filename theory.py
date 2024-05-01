import numpy as np
import math
from constants import Constants as Consts

class Theory:
    
    def __init__(self):
        self.consts = Consts()

    def Kn_density(self, T):
        '''
        Potassium number density
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
    
    def diamagnetic_FR(self, nu, l, T, B, nu_D1, nu_D2):
        """
        Calculations of diamagnetic Faraday rotation
        """
        Delta_Lambda_D1, Delta_nu_D1 = self.Zeeman_D1(B)
        Delta_Lambda_D2, Delta_nu_D2 = self.Zeeman_D2(B)

        print(f"Zeeman splitting of K D1 line with B={round(B*pow(10,4),3)} G is " , 
              str(Delta_Lambda_D1*pow(10,9)), " nm = ", str(Delta_nu_D1*pow(10,-9)), " GHz")
        print(f"Zeeman splitting of K D2 line with B={round(B*pow(10,4),3)} G is ",
              str(Delta_Lambda_D2*pow(10,9)), " nm = ", str(Delta_nu_D2*pow(10,-9)), " GHz")
        
        l = (7.5-0.159*2)*1e-2                                                                                                          # [m]

        term1 = ((7 / (nu - nu_D2)**2) + (4 / (nu - nu_D1)**2) - (2 / ((nu - nu_D2) * (nu - nu_D1)))) / (3 * self.consts.h)
        term2 = np.sign(B) * ((nu / (nu_D1 * (nu - nu_D1))) - (nu / (nu_D2 * (nu - nu_D2)))) / (self.consts.k_B * T)
        
        diamagnetic_theta = self.consts.alpha * self.consts.mu_B * self.Kn_density(T) * l  * B * (term1 + term2)

        return diamagnetic_theta
    
    def paramagnetic_FR(self, nu, l, T, P, nu_D1, nu_D2):
        """
        Calculations of paramagnetic Faraday rotation
        """
        l = (7.5-0.159*2)*1e-2                                                                                                          # [m]
        doppler_broad_D2 = self.doppler_broad(nu_D2, T)
        doppler_broad_D1 = self.doppler_broad(nu_D1, T)
        delta_nu_D2 = nu - nu_D2
        delta_nu_D1 = nu - nu_D1

        paramagnetic_thetea = (
            self.consts.alpha * self.Kn_density(T) * l * P * (
                (delta_nu_D2 / ((delta_nu_D2 - doppler_broad_D2) ** 2)) -
                (delta_nu_D1 / ((delta_nu_D1 - doppler_broad_D1) ** 2))
            )
        )                                                                                                                               # [rad]

        return paramagnetic_thetea
    
    def FR(self, nu, l, T, B, P, nu_D1, nu_D2):
        l = (7.5-0.159*2)*1e-2                                                                                                          # [m]
        
        theta = self.diamagnetic_FR(nu, l, T, B, nu_D1, nu_D2) + self.paramagnetic_FR(nu, l, T, P, nu_D1, nu_D2)

        return theta

    
    def Grad_theta(self, nu, l, T, B, nu_D1, nu_D2):
        """
        Wrong calculations
        """
        l = (7.5-0.159*2)*1e-2                                                                                                          # [m]

        # Gradient of Faraday rotation contributed from D1 line
        grad_D1_term1 = 8 * (Lambda**2/pow(Lambda-Lambda_D1,3) - Lambda/pow(Lambda-Lambda_D1,2)) / (3*self.Consts.c)
        grad_D1_term2 = -np.sign(B) * self.Consts.h / (self.Consts.k_B * (273.15+T) * pow(Lambda-Lambda_D1,2))
        grad_D1 = -self.Consts.alpha * self.Kn_density(T) * l * np.abs(B) * (Lambda_D1**2) * (grad_D1_term1 + grad_D1_term2) * 1e-6     # [microrad/pm]
        
        # Gradient of Faraday rotation contributed from D2 line
        grad_D2_term1 = 14 * (Lambda**2/pow(Lambda-Lambda_D2,3) - Lambda/pow(Lambda-Lambda_D2,2)) / (3*self.Consts.c)
        grad_D2_term2 = np.sign(B) * self.Consts.h / (self.Consts.k_B * (273.15+T) * pow(Lambda-Lambda_D2,2))
        grad_D2 = -self.Consts.alpha * self.Kn_density(T) * l * np.abs(B) * (Lambda_D1**2) * (grad_D2_term1 + grad_D2_term2) * 1e-6     # [microrad/pm]
        
        # Addition of gradient of Faraday rotation contributed from D1 & D2 line
        grad_theta = grad_D1 + grad_D2                                                                                                  # [microrad/pm]
        print('Grad_D1 is ' + str(grad_D1))
        print('Grad D2 is ' + str(grad_D2))
        
        return grad_theta