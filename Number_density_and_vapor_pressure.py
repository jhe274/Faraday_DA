import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

dir_path = os.path.join(
    'D:\\',  # Expands to your home directory
    'OneDrive', 
    'Files', 
    'Graduate_study', 
    'Research', 
    'Literature_Review',
    'Metal_Properties'
    )
fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

# Boltzmann Constant
kB = 1.3807 * ( 10 ** -23 )
mid_T = 450

# potassium number density
KT1 = np.arange(298.15, 336.5, 0.1)
KT2 = np.arange(336.5, mid_T-5, 0.1)
KT3 = np.arange(mid_T+5, 600, 0.1)

Kp1 = 10 ** ( 7.9667 - 4646 / KT1 ) # [Torr]
Kp2 = 10 ** ( 7.4077 - 4453 / KT2 ) # [Torr]
Kp3 = 10 ** ( 7.4077 - 4453 / KT3 ) # [Torr]
Kp_mid = 10 ** ( 7.4077 - 4453 / mid_T ) # [Torr]

rhoK1 = Kp1 / ((kB * KT1) * pow(10,6)) # [cm^-3]
rhoK2 = Kp2 / ((kB * KT2) * pow(10,6)) # [cm^-3]
rhoK3 = Kp3 / ((kB * KT3) * pow(10,6)) # [cm^-3]
rhoK_mid = Kp_mid / ((kB * mid_T) * pow(10,6)) # [cm^-3]

# rubidium number density
RbT1 = np.arange(298.15, 312.5, 0.1)
RbT2 = np.arange(312.5, mid_T-5, 0.1)
RbT3 = np.arange(mid_T+5, 600, 0.1)

Rbp1 = 10 ** ( 2.881 + 4.857 - 4215 / RbT1 ) # [Torr]
Rbp2 = 10 ** ( 2.881 + 4.312 - 4040 / RbT2 ) # [Torr]
Rbp3 = 10 ** ( 2.881 + 4.312 - 4040 / RbT3 ) # [Torr]
Rbp_mid = 10 ** ( 2.881 + 4.312 - 4040 / mid_T ) # [Torr]

rhoRb1 = Rbp1 / ((kB * RbT1) * pow(10,6))  # [cm^-3]
rhoRb2 = Rbp2 / ((kB * RbT2) * pow(10,6))  # [cm^-3]
rhoRb3 = Rbp3 / ((kB * RbT3) * pow(10,6))  # [cm^-3]
rhoRb_mid = Rbp_mid / ((kB * mid_T) * pow(10,6))  # [cm^-3]

# Caesium number density
CsT1 = np.arange(298.15, 301.7, 0.1)
CsT2 = np.arange(301.7, mid_T-5, 0.1)
CsT3 = np.arange(mid_T+5, 600, 0.1)

Csp1 = 10 ** ( -219.482 + 1088.676 / CsT1  - 0.08336185 * CsT1 + 94.88752 * np.log10(CsT1)) # [Torr]
Csp2 = 10 ** ( 8.22127 - 4006.048 / CsT2  - 0.00060194 * CsT2 - 0.19623 * np.log10(CsT2)) # [Torr]
Csp3 = 10 ** ( 8.22127 - 4006.048 / CsT3  - 0.00060194 * CsT3 - 0.19623 * np.log10(CsT3)) # [Torr]
Csp_mid = 10 ** ( 8.22127 - 4006.048 / mid_T  - 0.00060194 * mid_T - 0.19623 * np.log10(mid_T)) # [Torr]

rhoCs1 = Csp1 / ((kB * CsT1) * pow(10,6))  # [cm^-3]
rhoCs2 = Csp2 / ((kB * CsT2) * pow(10,6))  # [cm^-3]
rhoCs3 = Csp3 / ((kB * CsT3) * pow(10,6))  # [cm^-3]
rhoCs_mid = Csp_mid / ((kB * mid_T) * pow(10,6))  # [cm^-3]

# plot number density
ax.plot(KT1, rhoK1 * 133.322, color='blue', linewidth=2, label='Solid')
ax.plot(KT2, rhoK2 * 133.322, color='red', linewidth=2, label='Liquid')
ax.plot(KT3, rhoK3 * 133.322, color='red', linewidth=2)
ax.text(450, rhoK_mid * 133.322, 'K', fontsize=30, ha='center', va='center')
ax.plot(RbT1, rhoRb1 * 133.322, '--', color='blue', linewidth=2)
ax.plot(RbT2, rhoRb2 * 133.322, '--', color='red', linewidth=2)
ax.plot(RbT3, rhoRb3 * 133.322, '--', color='red', linewidth=2)
ax.text(450, rhoRb_mid * 133.322, 'Rb', fontsize=30, ha='center', va='center')
ax.plot(CsT1, rhoCs1 * 133.322, '-.', color='blue', linewidth=2)
ax.plot(CsT2, rhoCs2 * 133.322, '-.', color='red', linewidth=2)
ax.plot(CsT3, rhoCs3 * 133.322, '-.', color='red', linewidth=2)
ax.text(450, rhoCs_mid * 133.322, 'Cs', fontsize=30, ha='center', va='center')

# plot vapor pressure
# ax.plot(KT1, Kp1 , color='blue', linewidth=2, label='Solid')
# ax.plot(KT2, Kp2, color='red', linewidth=2, label='Liquid')
# ax.plot(KT3, Kp3, color='red', linewidth=2)
# ax.text(450, Kp_mid, 'K', fontsize=30, ha='center', va='center')
# ax.plot(RbT1, Rbp1, '--', color='blue', linewidth=2)
# ax.plot(RbT2, Rbp2, '--', color='red', linewidth=2)
# ax.plot(RbT3, Rbp3, '--', color='red', linewidth=2)
# ax.text(450, Rbp_mid, 'Rb', fontsize=30, ha='center', va='center')
# ax.plot(CsT1, Csp1, '-.', color='blue', linewidth=2)
# ax.plot(CsT2, Csp2, '-.', color='red', linewidth=2)
# ax.plot(CsT3, Csp3, '-.', color='red', linewidth=2)
# ax.text(450, Csp_mid, 'Cs', fontsize=30, ha='center', va='center')

ax.set_yscale("log", base=10)
ax.tick_params(axis='both', which='major', labelsize=30)
ax.legend(loc='best', fontsize=30)

# label the plot
ax.set_xlabel('Temperature, $T$ (K)', fontsize=30)
ax.set_ylabel('Number Density (cm$^{-3}$)', fontsize=30)
# ax.set_ylabel('Pressure (Torr)', fontsize=30)
# plt.title('Number density of potassium and rubidium')

plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(dir_path, 'Rb_K_Cs_Number_density.png'))
# plt.savefig(os.path.join(dir_path, 'Rb_K_Cs_Vapor_pressure.png'))
plt.show()