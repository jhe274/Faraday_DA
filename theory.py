import numpy as np
from constants import Constants as Consts

class Theory:
    
    def __init__(self):
        self.Consts = Consts()

    def Kn_density(self, T):
        '''
        Potassium number density
        '''
        if 25 < T <= 63.35:
            Kp = 10 ** (9.967 - 4646 / (273.15 + T))                                            # [Pa]
        elif 63.35 < T:
            Kp = 10 ** (9.408 - 4453 / (273.15 + T))                                            # [Pa]
        else:
            raise ValueError("Temperature T must be above 25°C")
        
        return Kp / (self.Consts.k_B * (273.15 + T))                                            # [m^-3]

    def Zeeman_ground(self, B):
        '''
        Calculations of Zeeman splitting 
        '''
        return self.Consts.g_G * (1/2) * self.Consts.mu_B * B

    def Zeeman_D1(self, B):
        """
        Zeeman splitting at potassium D1 line
        """
        E_D1 = self.Consts.g_D1 * (1/2) * self.Consts.mu_B * B
        Delta_E_D1 = -self.Zeeman_ground(B) - E_D1 - (self.Zeeman_ground(B) + E_D1)
        Delta_Lambda_D1 = - pow(self.Consts.Lambda_D1,2) * Delta_E_D1 / (self.Consts.h * self.Consts.c) # [m]
        Delta_nu_D1 = self.Consts.Nu_D1 - (self.Consts.c / (self.Consts.Lambda_D1 + Delta_Lambda_D1/2))*2 # [Hz]

        return Delta_Lambda_D1, Delta_nu_D1

    def Zeeman_D2(self, B):
        """
        Zeeman splitting at potassium D2 line
        """
        E_D2 = self.Consts.g_D2 * (3/2) * self.Consts.mu_B * B
        Delta_E_D2 = -self.Zeeman_ground(B) - E_D2 - (self.Zeeman_ground(B) + E_D2)
        Delta_Lambda_D2 = - pow(self.Consts.Lambda_D2,2) * Delta_E_D2 / (self.Consts.h * self.Consts.c) # [m]
        Delta_nu_D2 = self.Consts.Nu_D2 - (self.Consts.c / (self.Consts.Lambda_D2 + Delta_Lambda_D2/2)) * 2 # [Hz]

        return Delta_Lambda_D2, Delta_nu_D2
    
    def FR_theta(self, Lambda, l, B, T, Lambda_D1, Lambda_D2):
        """
        Calculations of polarization rotations
        """
        Delta_Lambda_D1, Delta_nu_D1 = self.Zeeman_D1(B)
        Delta_Lambda_D2, Delta_nu_D2 = self.Zeeman_D2(B)

        print(f"Zeeman splitting of K D1 line with B={round(B*pow(10,4),3)} G is " , 
              str(Delta_Lambda_D1*pow(10,9)), " nm = ", str(Delta_nu_D1*pow(10,-9)), " GHz")
        print(f"Zeeman splitting of K D2 line with B={round(B*pow(10,4),3)} G is ",
              str(Delta_Lambda_D2*pow(10,9)), " nm = ", str(Delta_nu_D2*pow(10,-9)), " GHz")
        
        l = (7.5-0.159*2)*1e-2                                                                  # [m]
        
        # Faraday rotation contributed from D1 line
        theta_D1_term1 = 4 * (Lambda**2) / (3 * self.Consts.c * pow(Lambda-Lambda_D1,2))
        theta_D1_term2 = -np.sign(B) * self.Consts.h / (self.Consts.k_B * (273.15+T) * (Lambda-Lambda_D1))
        theta_D1 = self.Consts.alpha * self.Kn_density(T) * l * np.abs(B) * (Lambda_D1**2) * (theta_D1_term1 + theta_D1_term2) # [rad]
        
        # Faraday rotation contributed from D2 line
        theta_D2_term1 = 7 * (Lambda**2) / (3 * self.Consts.c * pow(Lambda-Lambda_D2,2))
        theta_D2_term2 = np.sign(B) * self.Consts.h / (self.Consts.k_B * (273.15+T) * (Lambda-Lambda_D2))
        theta_D2 = self.Consts.alpha * self.Kn_density(T) * l * np.abs(B) * (Lambda_D2**2) * (theta_D2_term1 + theta_D2_term2) # [rad]
        
        # Addition of Faraday rotation contributed from D1 & D2 line
        theta = theta_D1 + theta_D2                                                             # [rad]
        
        return theta
    
    def Grad_theta(self, Lambda, B, T, Lambda_D1, Lambda_D2):
        """
        Calculations of change of polarization rotations w.r.t. wavelengths
        """
        l = (7.5-0.159*2)*1e-2                                                                  # [m]

        # Gradient of Faraday rotation contributed from D1 line
        grad_D1_term1 = 8 * (Lambda**2/pow(Lambda-Lambda_D1,3) - Lambda/pow(Lambda-Lambda_D1,2)) / (3*self.Consts.c)
        grad_D1_term2 = -np.sign(B) * self.Consts.h / (self.Consts.k_B * (273.15+T) * pow(Lambda-Lambda_D1,2))
        grad_D1 = -self.Consts.alpha * self.Kn_density(T) * l * np.abs(B) * (Lambda_D1**2) * (grad_D1_term1 + grad_D1_term2) * 1e-6 # [microrad/pm]
        
        # Gradient of Faraday rotation contributed from D2 line
        grad_D2_term1 = 14 * (Lambda**2/pow(Lambda-Lambda_D2,3) - Lambda/pow(Lambda-Lambda_D2,2)) / (3*self.Consts.c)
        grad_D2_term2 = np.sign(B) * self.Consts.h / (self.Consts.k_B * (273.15+T) * pow(Lambda-Lambda_D2,2))
        grad_D2 = -self.Consts.alpha * self.Kn_density(T) * l * np.abs(B) * (Lambda_D1**2) * (grad_D2_term1 + grad_D2_term2) * 1e-6 # [microrad/pm]
        
        # Addition of gradient of Faraday rotation contributed from D1 & D2 line
        grad_theta = grad_D1 + grad_D2                                                          # [microrad/pm]
        print('Grad_D1 is ' + str(grad_D1))
        print('Grad D2 is ' + str(grad_D2))
        
        return grad_theta