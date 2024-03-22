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
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()

    def XYplot(self, Bristol_t, Lambda, para, y_t, X, Y, run, n, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(2, 1, figsize=(25, 12))
        x, y1, y2 = [], [], []
        for i in range(run-1, run+3):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], X[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], X[i])
            Bristol_t[i], Lambda[i], y_t[i], Y[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], Y[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, X[i] = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], X[i][l_idx])
            Lambd, Y[i] = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], Y[i][l_idx])

            x.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                   # Frequency: [GHz]
            # y1.append(X[i] * 1e3)
            # y2.append(Y[i] * 1e3)
            y1.append(X[i] * 1e6)
            y2.append(Y[i] * 1e6)

        ax[0].plot(x[0], y1[0][1:], color='b', label=r'$\text{X}_\text{dc}$, $\theta_p$=89.5300°, $\theta_a$=46.2435°',
                   linestyle='-', linewidth=1, marker='^', markevery=200, markersize=10)
        ax[0].plot(x[0], y2[0][1:], color='b', label=r'$\text{Y}_\text{dc}$, $\theta_p$=89.5300°, $\theta_a$=46.2435°',
                    linestyle='-', linewidth=1, marker='x', markevery=200, markersize=10)
        ax[0].plot(x[1], y1[1][1:], color='r', label=r'$\text{X}_\text{dc}$, $\theta_p$=89.5070°, $\theta_a$=46.2365°',
                   linestyle='-', linewidth=1, marker='^', markevery=200, markersize=10)
        ax[0].plot(x[1], y2[1][1:], color='r', label=r'$\text{Y}_\text{dc}$, $\theta_p$=89.5070°, $\theta_a$=46.2365°',
                    linestyle='-', linewidth=1, marker='x', markevery=200, markersize=10)
        ax[0].tick_params(axis='both', labelsize=15)
        ax[0].set_xticks(np.arange(-5, 6, 1))
        ax[0].set_title(f'Vapor cell removed, run{run}-{run+1}', fontsize=20)
        ax[0].grid(True)
        ax[0].legend(loc='best', fontsize=15)

        # ax[0, 1].scatter(x[2], y1[2][1:], color='r', label=f'X{name}, L->H', marker='o', s=10)
        # ax[0, 1].scatter(x[2], y2[2][1:], color='r', label=f'Y{name}, L->H', marker='x', s=10)
        # ax[0, 1].scatter(x[3], y1[3][1:], color='b', label=f'X{name}, H->L', marker='o', s=10)
        # ax[0, 1].scatter(x[3], y2[3][1:], color='b', label=f'Y{name}, H->L', marker='x', s=10)
        # ax[0, 1].tick_params(axis='both', labelsize=15)
        # ax[0, 1].set_xticks(np.arange(-5, 6, 1))
        # ax[0, 1].set_title(r'Vapor cell inserted, $\lambda_{PEM}$=766.700 nm, $A$=2.391 rad', fontsize=20)
        # ax[0, 1].grid(True)
        # ax[0, 1].legend(loc='best', fontsize=15)

        ax[1].plot(x[2], y1[2][1:], color='b', label=r'$\text{X}_\text{dc}$, $\theta_p$=89.4835°, $\theta_a$=46.2290°', 
                   linestyle='-', linewidth=1, marker='^', markevery=200, markersize=10)
        ax[1].plot(x[2], y2[2][1:], color='b', label=r'$\text{Y}_\text{dc}$, $\theta_p$=89.4835°, $\theta_a$=46.2290°', 
                   linestyle='-', linewidth=1, marker='x', markevery=200, markersize=10)
        ax[1].plot(x[3], y1[3][1:], color='r', label=r'$\text{X}_\text{dc}$, $\theta_p$=89.4630°, $\theta_a$=46.2435°', 
                   linestyle='-', linewidth=1, marker='^', markevery=200, markersize=10)
        ax[1].plot(x[3], y2[3][1:], color='r', label=r'$\text{Y}_\text{dc}$, $\theta_p$=89.4630°, $\theta_a$=46.2435°', 
                   linestyle='-', linewidth=1, marker='x', markevery=200, markersize=10)
        ax[1].tick_params(axis='both', labelsize=15)
        ax[1].set_xticks(np.arange(-5, 6, 1))
        ax[1].set_title(f'Vapor cell removed, run{run+2}-{run+3}', fontsize=20)
        ax[1].grid(True)
        ax[1].legend(loc='best', fontsize=15)

        # ax[1, 1].scatter(x[6], y1[6][1:], color='r', label=f'X{name}, L->H', marker='o', s=10)
        # ax[1, 1].scatter(x[6], y2[6][1:], color='r', label=f'Y{name}, L->H', marker='x', s=10)
        # ax[1, 1].scatter(x[7], y1[7][1:], color='b', label=f'X{name}, H->L', marker='o', s=10)
        # ax[1, 1].scatter(x[7], y2[7][1:], color='b', label=f'Y{name}, H->L', marker='x', s=10)
        # ax[1, 1].tick_params(axis='both', labelsize=15)
        # ax[1, 1].set_xticks(np.arange(-5, 6, 1))
        # ax[1, 1].set_title(r'Vapor cell inserted, $\lambda_{PEM}$=766.690 nm, $A$=2.391 rad', fontsize=20)
        # ax[1, 1].grid(True)
        # ax[1, 1].legend(loc='best', fontsize=15)

        fig.text(0.5, 0.04, xlabel, ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle(title, fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'XY{name}_subplots_run{run}-{run+3}_@{date}.png'))
        plt.show()

    def Rplot(self, Bristol_t, Lambda, para, y_t, R, run, n, name, xlabel, ylabel, title):
        x, y = [], []
        fig, ax = plt.subplots(2, 2, figsize=(25, 12))  # 4 subplots, each with 1 column

        for i in range(run-1, run+3):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], R[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], R[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, V = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R[i][l_idx])
            x.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                         # Frequency: [GHz]
            y.append(V[1:] * 1e3)                                                                                                       # RMS Voltage: [mV]                            
            # y.append(V[1:] * 1e6)                                                                                                       # RMS Voltage: [microV]

        ax[0, 0].plot(x[0], y[0], label=r'$\theta_p$=89.5300°, $\theta_a$=46.2435°')
        # ax[0, 0].plot(x[1], y[1], label=r'$\theta_p$=89.5300°, $\theta_a$=46.2435°')
        ax[0, 0].tick_params(axis='both', labelsize=15)
        ax[0, 0].set_xticks(np.arange(-5, 6, 1))
        ax[0, 0].set_title(f'Vapor cell removed, run{run}', fontsize=20)
        ax[0, 0].grid(True)
        ax[0, 0].legend(loc='best', fontsize=15)

        ax[0, 1].plot(x[1], y[1], label=r'$\theta_p$=89.5070°, $\theta_a$=46.2365°')
        # ax[0, 1].plot(x[3], y[3], label=r'$\theta_p$=89.5070°, $\theta_a$=46.2365°')
        ax[0, 1].tick_params(axis='both', labelsize=15)
        ax[0, 1].set_xticks(np.arange(-5, 6, 1))
        ax[0, 1].set_title(f'Vapor cell removed, run{run+1}', fontsize=20)
        ax[0, 1].grid(True)
        ax[0, 1].legend(loc='best', fontsize=15)

        ax[1, 0].plot(x[2], y[2], label=r'$\theta_p$=89.4835°, $\theta_a$=46.2290°')
        # ax[1, 0].plot(x[5], y[5], label=r'$\theta_p$=89.4835°, $\theta_a$=46.2290°')
        ax[1, 0].tick_params(axis='both', labelsize=15)
        ax[1, 0].set_xticks(np.arange(-5, 6, 1))
        ax[1, 0].set_title(f'Vapor cell removed, run{run+2}', fontsize=20)
        ax[1, 0].grid(True)
        ax[1, 0].legend(loc='best', fontsize=15)

        ax[1, 1].plot(x[3], y[3], label=r'$\theta_p$=89.4630°, $\theta_a$=46.2435°')
        # ax[1, 1].plot(x[7], y[7], label=r'$\theta_p$=89.4630°, $\theta_a$=46.2435')
        ax[1, 1].tick_params(axis='both', labelsize=15)
        ax[1, 1].set_xticks(np.arange(-5, 6, 1))
        ax[1, 1].set_title(f'Vapor cell removed, run{run+3}', fontsize=20)
        ax[1, 1].grid(True)
        ax[1, 1].legend(loc='best', fontsize=15)

        fig.text(0.5, 0.04, xlabel, ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle(title, fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'R{name}_subplots_@{date}.png'))
        plt.show()

    def XY_vs_nu(self, lambda_path, lockins_path, name, run, n, B, P):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockins_path)

        if name == 'f':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, X1f, Y1f, run, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_{\omega}$ ($\mu$V)', r'$\text{XY}_{\omega}$ vs Frequency, run' + f'{run}-{run+3}' + 
                                  f', $B_z$={B} G, $P$={P} nW' + ' @'+ str(date))
        elif name == '2f':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, X2f, Y2f, run, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_{2\omega}$ ($\mu$V)', r'$\text{XY}_{2\omega}$ vs Frequency, run' + f'{run}-{run+3}' + 
                                  f', $B_z$={B} G, $P$={P} nW' + ' @'+ str(date))
        elif name == 'dc':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, Xdc, Ydc, run, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_\text{dc}$ (mV)', r'$\text{XY}_\text{dc}$ vs Frequency, run' + f'{run}-{run+3}' + 
                                  f', $B_z$={B} G, $P$={P} nW' + ' @'+ str(date))

    def R_vs_nu(self, lambda_path, lockins_path, name, run, n, B, P):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.FR_double_Kvapor(lockins_path)

        if name == 'f':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, R1f, run, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_{\omega}$ ($\mu$V)', r'$\text{R}_{\omega}$ vs Frequency, run' + f'{run}-{run+3}' + 
                                  f', $B_z$={B} G, $P$={P} nW' + ' @'+ str(date))
        elif name == '2f':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, R2f, run, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_{2\omega}$ ($\mu$V)', r'$\text{R}_{2\omega}$ vs Frequency, run' + f'{run}-{run+3}' + 
                                  f', $B_z$={B} G, $P$={P} nW' + ' @'+ str(date))
        elif name == 'dc':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, Rdc, run, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_\text{dc}$ (mV)', r'$\text{R}_\text{dc}$ vs Frequency, run' + f'{run}-{run+3}' + 
                                  f', $B_z$={B} G, $P$={P} nW' + ' @'+ str(date))

plotter = Plot()
date_input = '03-21-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.XY_vs_nu(Bristol_path, Lockins_path, 'f', 1, 5, 5.103, 865) 
# plotter.R_vs_nu(Bristol_path, Lockins_path, '2f', 1, 5, 5.103, 865) 