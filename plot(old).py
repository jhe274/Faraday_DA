import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

fig, ax = plt.subplots()
cwd = os.getcwd()
dir_path = os.path.join(cwd, 'Research', 'PhD Project', 'Faraday Rotation Measurements', 'K vapor cell')
Bristol = os.path.join(dir_path, 'Bristol data')
Lockins = os.path.join(dir_path, 'Lock-ins data')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()

    # Plot function of theoretical FR angle vs Wavelength
    def theta_theory_vs_lambda(self, l, B, T, Lambda_D1, Lambda_D2):
        for i in T:
            theta = Theory.FR_theta(l, B, i, Lambda_D1, Lambda_D2)
            Lambda = np.float64(np.arange(760, 777, 0.001) * pow(10,-9)) # [m]
            ax.plot(Lambda*pow(10,9), theta * pow(10,6), label=f'$T={i}°$C')
            # ax.plot(Lambda*pow(10,9), theta * pow(10,6) * 180 /np.pi, label=f'$T={i}°$C')

        ax.vlines(x=Consts.Lambda_D1*pow(10,9), ymin=0, ymax=325, colors='black', linestyles='--')
        # ax.vlines(x=Consts.Lambda_D2*pow(10,9), ymin=0, ymax=18500, colors='black', linestyles='--')
        # plt.annotate(r'$\lambda_{D1}$ =' f'{round(Consts.Lambda_D1*1e9,3)} nm', xy=(0.65,0.8), xycoords='axes fraction')
        plt.annotate(r'$\lambda_{D2}$ =' f'{round(Consts.Lambda_D2*1e9,3)} nm', xy=(0.25,0.8), xycoords='axes fraction')
        plt.xlabel(r'Wavelength (nm)')
        plt.ylabel(r'Polarization Rotation (microrad.)')
        # plt.ylabel(r'Polarization Rotation (microdeg.)')
        #ax.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.2f'))
        plt.xticks(np.arange(766.45,766.95,0.05))
        # plt.xticks(np.arange(round(CL_D2*1e9,3)-0.1, round(CL_D2*1e9,3)+0.1, 0.02))
        # plt.yticks(np.arange(0,18500,2000))
        plt.xlim(766.45,766.95)
        # plt.xlim(round(CL_D2*pow(10,9),3)-0.1, round(CL_D2*pow(10,9),3)+0.1)
        plt.ylim(0,325)
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'K Vapor Cell Faraday Rotation Angle vs Wavelength with $\\ell$={round(l*1e2,2)} cm, $B=${round(B*1e4,2)} G at $T={i}°$C Near D2 Line')
        plt.savefig(os.path.join(dir_path, 'Theory', 'Plots', 'K_vapor_cell_theta_vs_lambda_near_D_lines(rad)(test).png'))
        plt.show()

    # Plot function of measured temperature vs time of potassium vapor cell
    def temp_vs_time(self, temp_path, lock_in_path, i):
        Temp_t, T1, T2 = Read.TC300(temp_path)
        para, Lock_in_t, theta = Analyze.Triple_mod_theta(lock_in_path)

        diff = np.abs(Temp_t[i] - Lock_in_t[i][np.argmax(Lock_in_t[i])])
        end_t = np.argmin(diff)
        Temp_t[i] = Temp_t[i][:end_t+1]
        T1[i] = T1[i][:end_t+1]
        T2[i] = T2[i][:end_t+1]

        ax.plot(Temp_t[i]/60, T1[i], linewidth=1, label=r'Channel 1')
        ax.plot(Temp_t[i]/60, T2[i], linewidth=1, label=r'Channel 2')
        
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Temperature ($°$C)')
        ax.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.2f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Temperature vs Time of Potassium Vapor Cell @{date}')
        plt.savefig(os.path.join(Analysis, 'Plots', f'{date}', f'Temperature_vs_time_of_potassium_vapor_cell_{date}_run{i+1}.png'))
        plt.show()

    # Plot function of measured wavelength vs time
    def lambda_vs_time(self, lambda_path, lock_in_path, i, n):
        Bristol_t, Lambda = Read.Bristol(lambda_path)
        para, Lock_in_t, theta = Analyze.Triple_mod_theta(lock_in_path)
        
        Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i] = Analyze.filter_and_trim_data(Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i])
        idx = np.argmin(np.abs(Bristol_t[i] - (n+1)*1), axis=0)
        # print("The final time stamp of Bristol wavelength meter is " + str(Bristol_t[i][-1]))
        # print("The final time stamp of lock-in amplifiers is " + str(Lock_in_t[i][-1]))

        x = Bristol_t[i][idx:] # [s]
        y = Lambda[i][idx:] * 1e9 #[nm]
        y_ave = np.average(y)
        y_var = np.var(y, ddof=1) # unbiased sample variance
        y_std = np.std(y, ddof=1.5) # unbiased sample STD
        y_ste = y_std / np.sqrt(len(y))

        # Fit the measured wavelength with modified sine wave using curve fitting
        params, covariance = curve_fit(self.Analyze.modified_sine_wave, x, y, p0=[y_std*2, 10, 0, y_ave-y_std/2, 0])
        A_fit, f_fit, phi_fit, offset_fit, drift_slope_fit = params
        print('Delta lambda = ' + str(round(A_fit*2,6)) + ' [nm]')
        print('A = ' + str(A_fit))
        print('f = ' + str(f_fit))
        print('phi = ' + str(phi_fit))
        print('lambda_0 = ' + str(offset_fit))
        print('k = ' + str(drift_slope_fit))
        # plt.plot(x/60, self.Analyze.modified_sine_wave(x, A_fit, f_fit, phi_fit, offset_fit, drift_slope_fit), label='Fitted Sine Wave', color='red', linewidth=1)

        ax.scatter(x/60, y, s=.2, label=r'$\overline{{\lambda}}$={:.6f} nm, $\sigma$={:.6f} nm'.format(y_ave, y_std))
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Wavelength (nm)')
        ax.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.7f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Wavelength vs Time @{date}')
        plt.savefig(os.path.join(Analysis, 'Plots', f'{date}', f'Lambda_vs_time_{date}_run{i+1}.png'))
        plt.show()
    
    # Plot function of Vdc vs Time
    def Vdc_vs_time(self, lock_in_path, i):
        para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod = Read.Lock_ins(lock_in_path)

        # Calculate the values of Vdc
        Vdc = np.sqrt(Xdc[i] ** 2 + Ydc[i] ** 2) # [V]
        
        ax.scatter(Lock_in_t[i]/60, Vdc, s=0.5)
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'$V_{DC}$ (V)')
        plt.grid(True)
        # ax.legend(loc='best')
        plt.title(r'$V_{DC}$ vs Time @' + date)
        plt.savefig(os.path.join(Analysis, 'Plots', f'{date}', f'Vdc_vs_Time_{date}_run{i+1}.png'))
        plt.show()

    # Plot function of V2f vs Time
    def V2f_vs_time(self, lock_in_path, i):
        para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod = Read.Lock_ins(lock_in_path)

        # Calculate the values of V2f
        V2f = np.sqrt(X2f[i] ** 2 + Y2f[i] ** 2) # [V]
        
        ax.scatter(Lock_in_t[i]/60, V2f, s=0.5)
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'$V_{2f}$ (V)')
        plt.grid(True)
        # ax.legend(loc='best')
        plt.title(r'$V_{2f}$ vs Time @' + date)
        plt.savefig(os.path.join(Analysis, 'Plots', f'{date}', f'V2f_vs_Time_{date}_run{i+1}.png'))
        plt.show()
    
    # Plot function of Vmod vs Time
    def Vmod_vs_time(self, lock_in_path, i):
        para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod = Read.Lock_ins(lock_in_path)

        # Calculate the values of V2f
        Vmod = np.sqrt(Xmod[i] ** 2 + Ymod[i] ** 2) # [V]
        
        ax.scatter(Lock_in_t[i]/60, Vmod, s=2)
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'$V_{Mod}$ (V)')
        plt.grid(True)
        # ax.legend(loc='best')
        plt.title(r'$V_{Mod}$ vs Time @' + date)
        plt.savefig(os.path.join(Analysis, 'Plots', f'{date}', f'Vmod_vs_Time_{date}_run{i+1}.png'))
        plt.show()
    
    # Plot function of measured FR angle vs Time
    def theta_mea_vs_time(self, lambda_path, lock_in_path, i, n):
        Bristol_t, Lambda = Read.Bristol(lambda_path)
        para, Lock_in_t, theta = Analyze.Double_mod_theta(lock_in_path)
        # para, Lock_in_t, theta = Analyze.Triple_mod_theta(lock_in_path)

        # colors = plt.cm.viridis(np.linspace(0, 1, len(Lockins_path)))  # Generating a range of colors
        # color = colors[j]  # Selecting color for this plot

        Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i] = Analyze.filter_and_trim_data(Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i])
        Bristol_idx, lock_in_idx = Analyze.calculate_interval_and_indices(Bristol_t[i], Lock_in_t[i], para[i][2], n)
        
        x = Lock_in_t[i][lock_in_idx]
        y = theta[i][lock_in_idx] # [rad]
        y_ave = np.average(y)
        y_var = np.var(y, ddof=1) # unbiased sample variance
        y_std = np.std(y, ddof=1.5) # unbiased sample STD
        y_ste = y_std / np.sqrt(len(y))

        ax.scatter(x/60, y * pow(10,6), s=1,
                   label=r'$\overline{{\theta}}$={:.2f}$\pm${:.2f} $\mu$rad, $\sigma$={:.2f} $\mu$rad'.format(y_ave*1e6, y_ste*1e6, y_std*1e6))
            
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Polarization Rotation (microrad.)')
        plt.grid(True)
        # ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Time @{date}')
        plt.savefig(os.path.join(Analysis, 'Plots', f'{date}', f'Theta_vs_Time_filtered_{date}_run{i+1}.png'))
        plt.show()

    # Plot function of measured FR angle vs Wavelength
    def theta_mea_vs_lambda(self, lambda_path, lock_in_path, i, n):
        Bristol_t, Lambda = Read.Bristol(lambda_path)
        para, Lock_in_t, theta = Analyze.Double_mod_theta(lock_in_path)
        # para, Lock_in_t, theta = Analyze.Triple_mod_theta(lock_in_path)
        Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i] = Analyze.filter_and_trim_data(Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i])
        # itp_theta = Analyze.interp_theta(Bristol_t[i], Lock_in_t[i], theta[i])
        Bristol_idx, lock_in_idx = Analyze.calculate_interval_and_indices(Bristol_t[i], Lock_in_t[i], para[i][2], n)

        Lambda0 = np.linspace(766.696, 766.706, 10000) * 1e-9 # [m]
        Lambd = Lambda[i][Bristol_idx[:-50]]
        Thet = theta[i][lock_in_idx[:-50]]
        x = Lambd[np.nonzero(Lambd)] * 1e9 # [nm]
        y = Thet[np.nonzero(Thet)] * 1e3 # [millirad]
        x1 = x[:np.argmax(x)]
        x2 = x[np.argmax(x):]
        y1 = y[:np.argmax(x)]
        y2 = y[np.argmax(x):]
        valleys, _ = find_peaks(-y, prominence=.8, height=(-19, -15))
        ax.plot(x1, y1, color='black', label=r'forward scan')
        ax.plot(x2, y2, color='blue', label=r'backward scan')

        for j, k in enumerate(valleys):
            color = plt.cm.tab10(j)
            ax.scatter(x[k], y[k], color=color, s=20)
            ax.vlines(x=x[k], ymin=5, ymax=y[k], 
                      linestyle='--', color=color, label=r'$\lambda$={:.6f} nm, $\theta$={:.2f} millirad.'.format(x[k], y[k]), linewidth=2)
        
        plt.xlabel(r'Wavelength (nm)')
        plt.ylabel(r'Polarization Rotation (millirad.)')
        # plt.xticks(np.arange(766.696, 766.706, .001))
        # plt.yticks(np.arange(5,25,1))
        # plt.xlim(766.696, 766.706)
        # plt.ylim(5, 20)
        ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.6f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Wavelength, $B$=12.04 G, $P$=5 $\mu$W @{date}')
        plt.savefig(os.path.join(Analysis , 'Plots', f'{date}', f'Theta_vs_Lambda_{date}_run{i+1}.png'))
        plt.show()

    # Plot function of measured FR angle vs Frequency
    def theta_mea_vs_nu(self, lambda_path, lock_in_path, i, n):
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, Lock_in_t, theta = self.Analyze.Double_mod_theta(lock_in_path)
        # para, Lock_in_t, theta = Analyze.Triple_mod_theta(lock_in_path)
        Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
        Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i])
        l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], Lock_in_t[i], para[i][2], n)
        Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
        Thet = Thet[l_idx[1:-50]]

        print(np.where(np.isnan(Lambd[:-50])))
        # print(len(Lambd[:-50]), len(Thet[Thet != 0]))
        nu0 = self.Consts.c / self.Consts.Lambda_D2 * 1e-9 # [GHz]
        x = self.Consts.c / Lambd[:-50] * 1e-9 - nu0 # [MHz]
        y = Thet[Thet != 0] * 1e3 # [millirad]

        x1 = x[:np.argmin(x)]
        x2 = x[np.argmin(x):]
        y1 = y[:np.argmin(x)]
        y2 = y[np.argmin(x):]
        valleys, _ = find_peaks(-y, prominence=.5, height=(-45, -30))
        ax.plot(x1, y1, color='black', label=r'forward scan')
        ax.plot(x2, y2, color='blue', label=r'backward scan')

        for j, k in enumerate(valleys):
            color = plt.cm.tab10(j)
            ax.scatter(x[k], y[k], color=color, s=20)
            ax.vlines(x=x[k], ymin=5, ymax=y[k], 
                      linestyle='--', color=color, label=r'$\nu$={:.6f} GHz, $\theta$={:.2f} millirad.'.format(x[k], y[k]), linewidth=2)
        
        plt.xlabel(r'Frequency (GHz)')
        plt.ylabel(r'Polarization Rotation (millirad.)')
        plt.xticks(np.arange(self.Consts.c/766.705824 - nu0, self.Consts.c/766.696020 - nu0, .5))
        plt.yticks(np.arange(5,50,5))
        # plt.xlim(Consts.c/766.706, Consts.c/766.696)
        plt.ylim(5, 45)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Frequency, $B$=11.92 G, $P$=1.142 mW @{date}')
        plt.show()

    # Plot function of combined measured polarization rotations vs measured wavelengths
    def combined_theta_vs_lambda(self):
        Lambda = np.linspace(766.690, 766.712, 10000) * 1e-9 # [m]

        # Read data fom .csv file and convert them to numpy arrays
        data = Read.read_comb_theta_vs_lambda()
        Date, N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste = [
            np.array(arr) for arr in data
        ]
        # Sort by Lambda_ave and apply the same sorting to all arrays
        sort_idx = np.argsort(Lambda_ave)
        Date, N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste = [
            arr[sort_idx] for arr in
            (Date, N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste)
        ]

        x = [a*1e12 for a in Lambda_ave] # [pm]
        x_ste = [a*1e12 for a in Lambda_ste] # [pm]
        y = [(a / b)*1e-6 for a, b in zip(Theta_ave, np.abs(A_fit*2))] # [microrad/pm]
        y_std = [(a / b)*1e-6 for a, b in zip(Theta_std, np.abs(A_fit*2))] # [microrad/pm]
        y_ste = [(a / b)*1e-6 for a, b in zip(Theta_ste, np.abs(A_fit*2))] # [microrad/pm]
        ax.errorbar(x=x, xerr=x_ste, y=y, yerr=y_ste,
                fmt='o', markersize=3, elinewidth=2, capsize=3, color='red')  # mfc, mec

        params, covariance = curve_fit(Theory.Grad_theta, x, y, p0=None, 
                                    bounds=(( -9.7*1e-4, 26., Consts.CL_D1-0.1*1e-9, Consts.CL_D2-0.1*1e-9,), 
                                            ( -9.4*1e-4, 47., Consts.CL_D1+0.1*1e-9, Consts.CL_D2+0.1*1e-9,)),
                                    method='trf')
        B_fit, T_fit, D1_fit, D2_fit = params

        # Calculate the chi-squared value
        y_fit = Theory.Grad_theta(Lambda_ave, B_fit, T_fit, D1_fit, D2_fit)
        residuals = y - y_fit
        chi_squared = np.sum((residuals / y_std) ** 2)
        df = len(y) - len(params)
        reduced_chi_squared = chi_squared / df

        print('Fitted B-field is ' + str(B_fit))
        print('Fitted temperature is ' + str(T_fit))
        print('Fitted center of D1 line is ' + str(D1_fit))
        print('Fitted center of D2 line is ' + str(D2_fit))
        # print('Fitted alpha is ' + str(alpha_fit))
        # print('Fitted beta is ' + str(beta_fit))
        print('Fitted chi^2 is ' + str(chi_squared))
        grad_theta = Theory.Grad_theta(Lambda, B_fit, T_fit, D1_fit, D2_fit)
        ax.plot(Lambda*1e12, grad_theta, linestyle='--', label=f'$\chi^2$={chi_squared}')
        
        plt.xlabel(r'Wavelength (pm)')
        plt.ylabel(r'Polarization Rotation (microrad/pm)')
        ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        # plt.xticks(np.arange(766.65,766.75,0.01))
        # plt.yticks(np.linspace(0,2*1e6,10))
        plt.xlim(766690,766712)
        plt.ylim(-5*1e2,5*1e2)
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Wavelength')
        plt.savefig(os.path.join(Analysis, 'Plots', 'Measured_polarization_rotation_vs_wavelength_curve', 'Theta_vs_Lambda(test).png'))
        plt.show()

plotter = Plot()
date_input = input("Enter the date (MM-DD-YYYY): ")
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
# TC300_path = glob.glob(os.path.join(TC300, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.theta_mea_vs_nu(Bristol_path, Lockins_path, 8, 5)