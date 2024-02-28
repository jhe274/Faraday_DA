import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import scipy.special
from scipy.signal import find_peaks
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze
import numpy.polynomial.polynomial as poly

fig, ax = plt.subplots()
# dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
K_vapor = os.path.join(dir_path, 'K vapor cell')
Bristol = os.path.join(K_vapor, 'Bristol data')
Lockins = os.path.join(K_vapor, 'Lockins data')
DLCpro = os.path.join(K_vapor, 'TopticaDLCpro data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()

    def setup_plot(self, xlabel, ylabel, title, x):
        plt.xlabel(xlabel, fontsize=30)
        plt.ylabel(ylabel, fontsize=30)
        plt.title(title, fontsize=40)
        plt.grid(True)
        # plt.legend(loc='best')
        plt.xticks(np.linspace(x[1][0], x[1][-1], 10), fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlim(-0.8, 1.2)
        # plt.ylim(-1,1)
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True, useOffset=False))
        # plt.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))

    def myplot(self, x, y, label):
        ax.plot(x, y, label=label)

    def save_plot(self, name, i):
        plt.savefig(os.path.join(Plots, f'{date}', f'{name}_{date}_run{i+1}-{i+2}.png'))

    def process_and_plot(self, Bristol_t, Lambda, para, y_t, data, run, n, name, xlabel, ylabel, title):
        if np.all(para) == 0:
            # Wide Scan plot using DLC pro data
            x, y = [], []
            for i in range(run, run+2):
                Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
                Bristol_t[i], Lambda[i], y_t[i], data[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], y_t[i], data[i])
                l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], y_t[i], para, n)
                x.append(self.Consts.c / Lambda[i][b_idx] * 1e-9 - self.Consts.Nu39_D2 * 1e-9)              # [GHz]
                y.append(data[i][l_idx])
            
            # np.savetxt('data.txt', np.array([x[0],y[0]]).T)
            print(len(x[0]), len(y[0]))
            peaks, _ = find_peaks(y[0], prominence=.002, height=(0.57, .61))
            valleys, _ = find_peaks(-y[0], prominence=0.02, height=(-.57, -.5))
            K41_D2_A = (self.Consts.Nu41_D2_A-self.Consts.Nu39_D2) * 1e-9
            K41_D2_B = (self.Consts.Nu41_D2_B-self.Consts.Nu39_D2) * 1e-9
            K41_D2_C = (self.Consts.Nu41_D2_C-self.Consts.Nu39_D2) * 1e-9
            # Smooth the data using moving averages
            # window_size = 3
            # y_new = np.convolve(y[0], np.ones(window_size)/window_size, mode='valid')
            # x_new = x[0][window_size//2 : -window_size//2 + 1]                                         # Adjust x values to match smoothed data length
            derivative = np.gradient(y[0], x[0])

            fig, axs = plt.subplots(2, 1, figsize=(12, 12))
            axs[0].plot(x[0], y[0], color='red')
            for j,k in enumerate(peaks):
                axs[0].axvline(x=x[0][373], linestyle='--')
                axs[0].axvline(x=x[0][473], linestyle='--')
            for j,k in enumerate(valleys):
                axs[0].axvline(x=x[0][k], linestyle='--')
            axs[0].axvline(x=K41_D2_A, color='orange', linestyle='--')
            axs[0].axvline(x=K41_D2_B, color='orange', linestyle='--')
            axs[0].axvline(x=K41_D2_C, color='orange', linestyle='--')
            axs[0].grid(True)
            axs[0].set_xlim(-0.8, 1.2)
            # axs[0].xaxis.set_major_formatter(ScalarFormatter(useMathText=True, useOffset=False))
            axs[0].tick_params(axis='x', labelsize=15)
            axs[0].tick_params(axis='y', labelsize=15)
            axs[0].set_ylabel('Intensity (au)', fontsize=20)
            axs[0].set_xlabel('Frequency (GHz)', fontsize=20)
            
            axs[1].plot(x[0], derivative, color='blue')
            for j,k in enumerate(peaks):
                axs[1].axvline(x=x[0][373], linestyle='--')
                axs[1].axvline(x=x[0][473], linestyle='--')
            for j,k in enumerate(valleys):
                axs[1].axvline(x=x[0][k], linestyle='--')
            axs[1].grid(True)
            axs[1].set_xlim(-0.8, 1.2)
            axs[1].set_ylim(-2.2,1.5)
            axs[1].tick_params(axis='x', labelsize=15)
            axs[1].tick_params(axis='y', labelsize=15)
            axs[1].set_ylabel('Intensity (au)', fontsize=20)
            axs[1].set_xlabel('Frequency (GHz)', fontsize=20)

        else:
            # RMS plot using lock-ins data
            for i in range(run, run+2):
                Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
                Bristol_t[i], Lambda[i], y_t[i], data[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], y_t[i], data[i])
                l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
                Lambd, y = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], data[i][l_idx])
                x = self.Consts.c / Lambd * 1e-9 - self.Consts.Nu39_D2 * 1e-9                         # [GHz]

                # Baseline fit
                # coefs = poly.polyfit(x, y[1:], 2)
                # ffit = poly.polyval(x, coefs)

                self.scatter_plot(x, y[1:], f'{r"H->L" if i == run else r"L->H"} Wide Scan')
                # plt.plot(x, ffit, ':', label='Baseline Fit')

        # plt.axvline(x=0, color = 'red', linestyle='--')
        # self.setup_plot(xlabel, ylabel, title, x)
        plt.tight_layout()  # Ensure tight layout to prevent overlapping
        self.save_plot(name, i)
        # plt.show()
        plt.close()

    def y_vs_nu(self, lambda_path, lockins_path, y_path, name, run, n, gain, B, power):
        """
        Plot function of Lock-in RMS values vs Wavelength
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, Rmod, R2f, Rdc, _ = self.Analyze.Triple_modu_theta(lockins_path)
        x, y, Y, DLCpro_t = self.Read.DLCpro_WideScan(y_path)
        if name == 'SAS_KVC':
            self.process_and_plot(Bristol_t, Lambda, 0, DLCpro_t, y, run-1, n, name, 'Frequency (GHz)',
                                #   'Intensity (au)', r'Doppler-free SAS of K vapor cell @$T=36°$C')
                                  'Intensity (au)', r'Doppler-free SAS of K vapor cell with blanced photodetector @$T=36°$C')
                                #   r'$V$ (V)', r'Saturated absorption spectrum of K vapor cell, run' + f'{run}-{run+1}' + ' @'+ str(date))
        elif name == 'Rmod':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, Rmod, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{Mod}$ (V)', r'$R_{Mod}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', PD gain={gain}, $|B|$={B} G, $P$={power} $\mu$W' + ' @'+ str(date))
        elif name == 'R2f':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, R2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{2f}$ (V)', r'$R_{2f}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', PD gain={gain}, $|B|$={B} G, $P$={power} $\mu$W' + ' @'+ str(date))
        elif name == 'Rdc':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, Rdc, run-1, n, name, 'Frequency GHz)',
                                  r'$R_{DC}$ (V)', r'$R_{DC}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', PD gain={gain}, $|B|$={B} G, $P$={power} $\mu$W' + ' @'+ str(date))
    
    def frequency_detuning(self, lambda_path, lockin_path, run, n, gain, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, R2f, Rdc, theta = self.Analyze.Double_modu_theta(lockin_path)
        xmin, xmax = float('-inf'), float('inf')
        fit_params = []
        for i in range(run-1, run+1):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])

            # Fit
            coefs = poly.polyfit(lockins_t[i], Rdc[i], 2)
            Rdc_fit = poly.polyval(lockins_t[i], coefs)
            theta_fit = R2f[i] / (np.pi * scipy.special.jv(2,2.405) * Rdc_fit)
            Bristol_t[i], Lambda[i], lockins_t[i], theta_fit = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta_fit)
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd_fit, Thet_fit = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta_fit[l_idx])
            x_fit = self.Consts.c / Lambd_fit * 1e-9 - self.Consts.Nu39_D2 * 1e-9                     # [GHz]
            y_fit = Thet_fit[1:] * 1e3                                                              # [millirad]
            

            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x = self.Consts.c / Lambd * 1e-9 - self.Consts.Nu39_D2 * 1e-9                             # [GHz]
            y = Thet[1:] * 1e3                                                                      # [millirad]
            
            # valleys, _ = find_peaks(-y, prominence=.5, height=(-28, -23))
            # color = 'black' if i == run-1 else 'blue'
            # TC = ['5ms', '10ms', '20ms', '50ms', '100ms']
            ax.plot(x, y, label=f'{r"H->L" if i == run-1 else r"L->H"} Wide Scan')
            # plt.plot(x_fit, y_fit, ':', label=f'{r"H->L" if i == 0 else r"L->H"} Wide Scan')
            # ax.plot(x, y, label=f'TC = {TC[i-5]}')

            xmin = max(xmin, min(x))
            xmax = min(xmax, max(x))

        # lambda_theo = np.linspace(766.695*1e-9, 766.705*1e-9, 2000)
        # theta_theo = self.Theory.FR_theta(lambda_theo, 0.0718, 5*1e-4, 35, self.Consts.Lambda39_D1, self.Consts.Lambda39_D2)
        # x_theo = self.Consts.c / lambda_theo * 1e-9 - self.Consts.Nu39_D2 
        # ax.plot(x_theo, theta_theo * 1e3, '--', color='red', label='Theory')

        # for j, k in enumerate(valleys):
        #     color = plt.cm.tab10(j)
        #     ax.scatter(x[k], y[k], color=color, s=20)
            # ax.vlines(x=x[k], ymin=16, ymax=y[k], 
                    #   linestyle='--', color=color, label=r'$\nu$={:.6f} GHz, $\theta$={:.2f} millirad.'.format(x[k], y[k]), linewidth=2)
        
        plt.xlabel(r'Frequency (GHz)')
        plt.ylabel(r'Polarization Rotation (millirad.)')
        plt.xticks(np.arange(np.floor(xmin), np.ceil(xmax)+1, 0.5))
        # plt.yticks(np.arange(5,50,5))
        # plt.ylim(0, 22)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Frequency, run{run}-{run+1}, PD gain={gain}, $|B|$={B} G, $P$={power} $\mu$W @{date}')
        plt.show()

plotter = Plot()
# date_input = input("Enter the date (MM-DD-YYYY): ")
date_input = '02-13-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
DLCpro_path = glob.glob(os.path.join(DLCpro, date, '*.csv'))
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
# plotter.frequency_detuning(Bristol_path, Lockins_path, 1, 5, 6, 5.17, 2.28)
plotter.y_vs_nu(Bristol_path, Lockins_path, DLCpro_path, 'SAS_KVC', 1, 5, 6, 5.16, 2.28)