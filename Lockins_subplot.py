import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
# dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
K_vapor = os.path.join(dir_path, 'K vapor cell')
Bristol = os.path.join(K_vapor, 'Bristol data')
Lockins = os.path.join(K_vapor, 'Lockins data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()

    def process_and_plot(self, Bristol_t, Lambda, para, y_t, R, run, n, name, xlabel, ylabel, title):
        x, y = [], []
        fig, ax = plt.subplots(2, 2, figsize=(25, 12))  # 4 subplots, each with 1 column

        for i in range(run, run+8):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], R[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], y_t[i], R[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, V = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R[i][l_idx])
            x.append(self.Consts.c / Lambd * 1e-9 - self.Consts.Nu39_D2 * 1e-9)                                                         # Frequency: [GHz]
            # y.append(V[1:] * 1e3)                                                                                                       # RMS Voltage: [mV]                            
            y.append(V[1:] * 1e6)                                                                                                       # RMS Voltage: [microV]

        ax[0, 0].plot(x[0], y[0], label=f'{"L->H"} Wide Scan')
        ax[0, 0].plot(x[1], y[1], label=f'{"H->L"} Wide Scan')
        ax[0, 0].tick_params(axis='x', labelsize=15)
        ax[0, 0].tick_params(axis='y', labelsize=15)
        ax[0, 0].set_title(r'Vapor cell removed, $P$=520 nW, $\lambda_{PEM}$=766.701 nm', fontsize=20)
        ax[0, 0].grid(True)
        ax[0, 0].legend(loc='best')

        ax[0, 1].plot(x[2], y[2], label=f'{"H->L"} Wide Scan')
        ax[0, 1].plot(x[3], y[3], label=f'{"L->H"} Wide Scan')
        ax[0, 1].tick_params(axis='x', labelsize=15)
        ax[0, 1].tick_params(axis='y', labelsize=15)
        ax[0, 1].set_title(r'Vapor cell inserted, $P$=570 nW, $\lambda_{PEM}$=766.701 nm', fontsize=20)
        ax[0, 1].grid(True)
        ax[0, 1].legend(loc='best')

        ax[1, 0].plot(x[4], y[4], label=f'{"L->H"} Wide Scan')
        ax[1, 0].plot(x[5], y[5], label=f'{"H->L"} Wide Scan')
        ax[1, 0].tick_params(axis='x', labelsize=15)
        ax[1, 0].tick_params(axis='y', labelsize=15)
        ax[1, 0].set_title(r'Vapor cell removed, $P$=550 nW, $\lambda_{PEM}$=766.691 nm', fontsize=20)
        ax[1, 0].grid(True)
        ax[1, 0].legend(loc='best')

        ax[1, 1].plot(x[6], y[6], label=f'{"L->H"} Wide Scan')
        ax[1, 1].plot(x[7], y[7], label=f'{"H->L"} Wide Scan')
        ax[1, 1].tick_params(axis='x', labelsize=15)
        ax[1, 1].tick_params(axis='y', labelsize=15)
        ax[1, 1].set_title(r'Vapor cell inserted, $P$=570 nW, $\lambda_{PEM}$=766.691 nm', fontsize=20)
        ax[1, 1].grid(True)
        ax[1, 1].legend(loc='best')

        fig.text(0.5, 0.04, xlabel, ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle(title, fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'{name}_subplots_@{date}.png'))
        plt.show()


    def y_vs_nu(self, lambda_path, lockins_path, name, run, n, B, T):
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.Analyze.Double_modu_theta(lockins_path)

        if name == 'R1f':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, R1f, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{1f}$ ($\mu$V)', r'$R_{1f}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $T$={T}°C' + ' @'+ str(date))
        elif name == 'R2f':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, R2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{2f}$ ($\mu$V)', r'$R_{2f}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $T$={T}°C' + ' @'+ str(date))
        elif name == 'Rdc':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, Rdc, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{dc}$ (mV)', r'$R_{dc}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $T$={T}°C' + ' @'+ str(date))

plotter = Plot()
date_input = '03-01-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.y_vs_nu(Bristol_path, Lockins_path, 'R1f', 1, 5, 5.103, 25) 