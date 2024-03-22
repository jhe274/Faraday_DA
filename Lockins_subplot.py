import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

# dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
K_vapor = os.path.join(dir_path, 'K vapor cell')
Bristol = os.path.join(K_vapor, 'Bristol data')
Lockins = os.path.join(K_vapor, 'Lockins data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()

    def XYplot(self, Bristol_t, Lambda, para, y_t, X, Y, run, n, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(2, 2, figsize=(25, 12))
        x, y1, y2 = [], [], []
        for i in range(run, run+8):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], X[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], X[i])
            Bristol_t[i], Lambda[i], y_t[i], Y[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], Y[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, X[i] = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], X[i][l_idx])
            Lambd, Y[i] = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], Y[i][l_idx])

            x.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                   # Frequency: [GHz]
            y1.append(X[i] * 1e3)
            y2.append(Y[i] * 1e3)
            # y1.append(X[i] * 1e6)
            # y2.append(Y[i] * 1e6)

        ax[0, 0].scatter(x[0], y1[0][1:], color='r', label=f'X{name}, L->H', marker='o', s=10)
        ax[0, 0].scatter(x[0], y2[0][1:], color='r', label=f'Y{name}, L->H', marker='x', s=10)
        ax[0, 0].scatter(x[1], y1[1][1:], color='b', label=f'X{name}, H->L', marker='o', s=10)
        ax[0, 0].scatter(x[1], y2[1][1:], color='b', label=f'Y{name}, H->L', marker='x', s=10)
        ax[0, 0].tick_params(axis='both', labelsize=15)
        ax[0, 0].set_xticks(np.arange(-5, 6, 1))
        ax[0, 0].set_title(r'Vapor cell removed, $\lambda_{PEM}$=766.700 nm, $A$=2.391 rad', fontsize=20)
        ax[0, 0].grid(True)
        ax[0, 0].legend(loc='best', fontsize=15)

        ax[0, 1].scatter(x[2], y1[2][1:], color='r', label=f'X{name}, L->H', marker='o', s=10)
        ax[0, 1].scatter(x[2], y2[2][1:], color='r', label=f'Y{name}, L->H', marker='x', s=10)
        ax[0, 1].scatter(x[3], y1[3][1:], color='b', label=f'X{name}, H->L', marker='o', s=10)
        ax[0, 1].scatter(x[3], y2[3][1:], color='b', label=f'Y{name}, H->L', marker='x', s=10)
        ax[0, 1].tick_params(axis='both', labelsize=15)
        ax[0, 1].set_xticks(np.arange(-5, 6, 1))
        ax[0, 1].set_title(r'Vapor cell inserted, $\lambda_{PEM}$=766.700 nm, $A$=2.391 rad', fontsize=20)
        ax[0, 1].grid(True)
        ax[0, 1].legend(loc='best', fontsize=15)

        ax[1, 0].scatter(x[4], y1[4][1:], color='r', label=f'X{name}, L->H', marker='o', s=10)
        ax[1, 0].scatter(x[4], y2[4][1:], color='r', label=f'Y{name}, L->H', marker='x', s=10)
        ax[1, 0].scatter(x[5], y1[5][1:], color='b', label=f'X{name}, H->L', marker='o', s=10)
        ax[1, 0].scatter(x[5], y2[5][1:], color='b', label=f'Y{name}, H->L', marker='x', s=10)
        ax[1, 0].tick_params(axis='both', labelsize=15)
        ax[1, 0].set_xticks(np.arange(-5, 6, 1))
        ax[1, 0].set_title(r'Vapor cell removed, $\lambda_{PEM}$=766.690 nm, $A$=2.391 rad', fontsize=20)
        ax[1, 0].grid(True)
        ax[1, 0].legend(loc='best', fontsize=15)

        ax[1, 1].scatter(x[6], y1[6][1:], color='r', label=f'X{name}, L->H', marker='o', s=10)
        ax[1, 1].scatter(x[6], y2[6][1:], color='r', label=f'Y{name}, L->H', marker='x', s=10)
        ax[1, 1].scatter(x[7], y1[7][1:], color='b', label=f'X{name}, H->L', marker='o', s=10)
        ax[1, 1].scatter(x[7], y2[7][1:], color='b', label=f'Y{name}, H->L', marker='x', s=10)
        ax[1, 1].tick_params(axis='both', labelsize=15)
        ax[1, 1].set_xticks(np.arange(-5, 6, 1))
        ax[1, 1].set_title(r'Vapor cell inserted, $\lambda_{PEM}$=766.690 nm, $A$=2.391 rad', fontsize=20)
        ax[1, 1].grid(True)
        ax[1, 1].legend(loc='best', fontsize=15)

        fig.text(0.5, 0.04, xlabel, ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle(title, fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'{name}_subplots_@{date}.png'))
        plt.show()

    def Rplot(self, Bristol_t, Lambda, para, y_t, R, run, n, name, xlabel, ylabel, title):
        x, y = [], []
        fig, ax = plt.subplots(2, 2, figsize=(25, 12))  # 4 subplots, each with 1 column

        for i in range(run, run+8):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], R[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], R[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, V = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R[i][l_idx])
            x.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                         # Frequency: [GHz]
            y.append(V[1:] * 1e3)                                                                                                       # RMS Voltage: [mV]                            
            # y.append(V[1:] * 1e6)                                                                                                       # RMS Voltage: [microV]

        ax[0, 0].plot(x[0], y[0], label=f'{"L->H"} Wide Scan')
        ax[0, 0].plot(x[1], y[1], label=f'{"H->L"} Wide Scan')
        ax[0, 0].tick_params(axis='both', labelsize=15)
        ax[0, 0].set_xticks(np.arange(-5, 6, 1))
        ax[0, 0].set_title(r'Vapor cell removed, $P$=875 nW, $\lambda_{PEM}$=766.700 nm, $A$=2.391 rad', fontsize=20)
        ax[0, 0].grid(True)
        ax[0, 0].legend(loc='best', fontsize=15)

        ax[0, 1].plot(x[2], y[2], label=f'{"L->H"} Wide Scan')
        ax[0, 1].plot(x[3], y[3], label=f'{"H->L"} Wide Scan')
        ax[0, 1].tick_params(axis='both', labelsize=15)
        ax[0, 1].set_xticks(np.arange(-5, 6, 1))
        ax[0, 1].set_title(r'Vapor cell inserted, $P$=870 nW, $\lambda_{PEM}$=766.700 nm, $A$=2.391 rad', fontsize=20)
        ax[0, 1].grid(True)
        ax[0, 1].legend(loc='best', fontsize=15)

        ax[1, 0].plot(x[4], y[4], label=f'{"L->H"} Wide Scan')
        ax[1, 0].plot(x[5], y[5], label=f'{"H->L"} Wide Scan')
        ax[1, 0].tick_params(axis='both', labelsize=15)
        ax[1, 0].set_xticks(np.arange(-5, 6, 1))
        ax[1, 0].set_title(r'Vapor cell removed, $P$=876 nW, $\lambda_{PEM}$=766.690 nm, $A$=2.391 rad', fontsize=20)
        ax[1, 0].grid(True)
        ax[1, 0].legend(loc='best')

        ax[1, 1].plot(x[6], y[6], label=f'{"L->H"} Wide Scan')
        ax[1, 1].plot(x[7], y[7], label=f'{"H->L"} Wide Scan')
        ax[1, 1].tick_params(axis='both', labelsize=15)
        ax[1, 1].set_xticks(np.arange(-5, 6, 1))
        ax[1, 1].set_title(r'Vapor cell inserted, $P$=873 nW, $\lambda_{PEM}$=766.690 nm, $A$=2.391 rad', fontsize=20)
        ax[1, 1].grid(True)
        ax[1, 1].legend(loc='best')

        fig.text(0.5, 0.04, xlabel, ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle(title, fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'{name}_subplots_@{date}.png'))
        plt.show()

    def XY_vs_nu(self, lambda_path, lockins_path, name, run, n, B, power):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockins_path)

        if name == '1f':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, X1f, Y1f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_\text{1f}$ ($\mu$V)', r'$\text{XY}_\text{1f}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $P$={power} nW' + ' @'+ str(date))
        elif name == '2f':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, X2f, Y2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_\text{2f}$ ($\mu$V)', r'$\text{XY}_\text{2f}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $P$={power} nW' + ' @'+ str(date))
        elif name == 'dc':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, Xdc, Ydc, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_\text{dc}$ (mV)', r'$\text{XY}_\text{dc}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $P$={power} nW' + ' @'+ str(date))

    def R_vs_nu(self, lambda_path, lockins_path, name, run, n, B, T):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.FR_double_Kvapor(lockins_path)

        if name == '1f':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, R1f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_\text{1f}$ ($\mu$V)', r'$\text{R}_\text{1f}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $T$={T}°C' + ' @'+ str(date))
        elif name == '2f':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, R2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_\text{2f}$ ($\mu$V)', r'$\text{R}_\text{2f}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $T$={T}°C' + ' @'+ str(date))
        elif name == 'dc':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, Rdc, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_\text{dc}$ (mV)', r'$\text{R}_\text{dc}$ vs Frequency, run' + f'{run}-{run+7}' + 
                                  f', $B_z$={B} G, $T$={T}°C' + ' @'+ str(date))

plotter = Plot()
date_input = '03-19-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.XY_vs_nu(Bristol_path, Lockins_path, '1f', 1, 5, 5.103, 870) 
# plotter.R_vs_nu(Bristol_path, Lockins_path, '2f', 1, 5, 5.103, 25) 