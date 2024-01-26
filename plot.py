import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.signal import find_peaks
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

fig, ax = plt.subplots()
dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
K_vapor = os.path.join(dir_path, 'K vapor cell')
# dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements', 'K vapor cell')
Bristol = os.path.join(K_vapor, 'Bristol data')
Lockins = os.path.join(K_vapor, 'Lockins data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()

    def setup_plot(self, xlabel, ylabel, title, x):
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.xticks(np.linspace(x[0], x[-1], 10))
        # plt.xlim(x[0], x[-1])
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True, useOffset=False))
        # plt.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))

    def scatter_plot(self, x, y, label, s=2):
        ax.scatter(x, y, s=s, label=label)

    def save_plot(self, name, i):
        plt.savefig(os.path.join(Plots, f'{date}', f'R{name}_{date}_run{i+1}.png'))

    def process_and_plot(self, Bristol_t, Lambda, para, lockins_t, data, i, n, name, xlabel, ylabel, title):
        Bristol_t[i], Lambda[i], lockins_t[i], data[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], data[i])
        l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
        Lambd, R = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], data[i][l_idx])

        self.scatter_plot(Lambd * 1e9, R[1:], name)
        self.setup_plot(xlabel, ylabel, title, Lambd * 1e9)
        self.save_plot(name, i)
        plt.show()

    def R_vs_wvl(self, lambda_path, lockin_path, name, i, n):
        """
        Plot function of Lock-in RMS values vs Wavelength
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, Rmod, R2f, Rdc, _ = self.Analyze.Triple_modu_theta(lockin_path)

        if name == 'mod':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, Rmod, i-1, n, name, 'Wavelength (nm)',
                                  r'$R_{Mod}$', r'$R_{Mod}$ vs Wavelength, run' + str(i) + ' @'+ str(date))
        elif name == '2f':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, R2f, i-1, n, name, 'Wavelength (nm)',
                                  r'$R_{2f}$', r'$R_{2f}$ vs Wavelength, run' + str(i) + ' @'+ str(date))
        elif name == 'dc':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, Rdc, i-1, n, name, 'Wavelength (nm)',
                                  r'$R_{DC}$', r'$R_{DC}$ vs Wavelength, run' + str(i) + ' @'+ str(date))

    def frequency_detuning(self, lambda_path, lockin_path, n):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, R2f, Rdc, theta = self.Analyze.Double_modu_theta(lockin_path)
        xmin, xmax = float('-inf'), float('inf')
        for i in range(6,8):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            nu0 = self.Consts.c / self.Consts.Lambda_D2 * 1e-9                              # [GHz]
            x = self.Consts.c / Lambd * 1e-9 - nu0                                          # [GHz]
            y = Thet[1:] * 1e3                                                              # [millirad]

            valleys, _ = find_peaks(-y, prominence=.5, height=(-28, -23))
            color = 'black' if i == 6 else 'blue'
            ax.plot(x, y, color=color, label=f'{r"L->H" if i == 6 else r"H->L"} Wide Scan')

            xmin = max(xmin, min(x))
            xmax = min(xmax, max(x))

        lambda_theo = np.linspace(766.695*1e-9, 766.705*1e-9, 2000)
        theta_theo = self.Theory.FR_theta(lambda_theo, 0.0718, 5*1e-4, 35, self.Consts.Lambda_D1, self.Consts.Lambda_D2)
        x_theo = self.Consts.c / lambda_theo * 1e-9 - nu0
        ax.plot(x_theo, theta_theo * 1e3, '--', color='red', label='Theory')

        for j, k in enumerate(valleys):
            color = plt.cm.tab10(j)
            ax.scatter(x[k], y[k], color=color, s=20)
            ax.vlines(x=x[k], ymin=16, ymax=y[k], 
                      linestyle='--', color=color, label=r'$\nu$={:.6f} GHz, $\theta$={:.2f} millirad.'.format(x[k], y[k]), linewidth=2)
        
        plt.xlabel(r'Frequency (GHz)')
        plt.ylabel(r'Polarization Rotation (millirad.)')
        plt.xticks(np.arange(np.floor(xmin), np.ceil(xmax), 0.5))
        # plt.yticks(np.arange(5,50,5))
        plt.ylim(0, 22)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Frequency, $B$=5.15 G, $P$=2.64 $\mu$W @{date}')
        plt.show()

plotter = Plot()
# date_input = input("Enter the date (MM-DD-YYYY): ")
date_input = '01-25-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
# plotter.frequency_detuning(Bristol_path, Lockins_path, 5)
plotter.R_vs_wvl(Bristol_path, Lockins_path, 'dc', 6, 5)