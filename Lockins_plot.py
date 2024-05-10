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
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in range(run, run+2):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], X[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], X[i])
            Bristol_t[i], Lambda[i], y_t[i], Y[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], Y[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, X[i] = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], X[i][l_idx])
            Lambd, Y[i] = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], Y[i][l_idx])
            x = self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9                                                   # Frequency: [GHz]
            colors = 'y' if i == run else ('r' if i == run+1 else 'b')

            if name == 'f':
                scale_factor = 1e6  # RMS Voltage: [microV]
            elif name == '2f':
                scale_factor = 1e6  # RMS Voltage: [microV]
            else:
                scale_factor = 1e3  # RMS Voltage: [mV]

            # Apply the scaling factor to X and Y arrays
            X[i] = X[i] * scale_factor
            Y[i] = Y[i] * scale_factor

            # Plot X and Y with the appropriate labels
            label_x = (r'$\text{X}_\text{f}$' if name == 'f' else 
                    r'$\text{X}_\text{2f}$' if name == '2f' else 
                    r'$\text{X}_\text{dc}$') + (', air' if i == run else ', vapor cell')
            label_y = (r'$\text{Y}_\text{f}$' if name == 'f' else 
                    r'$\text{Y}_\text{2f}$' if name == '2f' else 
                    r'$\text{Y}_\text{dc}$') + (', air' if i == run else ', vapor cell')

            ax.plot(x, X[i][1:], color=colors, label=label_x, linestyle='-', linewidth=1, marker='^', markevery=200, markersize=10)
            ax.plot(x, Y[i][1:], color=colors, label=label_y, linestyle='-', linewidth=1, marker='x', markevery=200, markersize=10)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'XY{name}_{date}_run{i}-{i+1}.png'))
        plt.show()

    def Rplot(self, Bristol_t, Lambda, para, y_t, R, run, n, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        for i in range(run, run+2):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], R[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], R[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, y = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R[i][l_idx])
            x = self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9                                                   # Frequency: [GHz]
            # y = y * 1e3                                                                                                     # RMS Voltage: [mV]                            
            y = y * 1e6                                                                                                     # RMS Voltage: [microV]

            ax.scatter(x, y[1:], label="cell inserted" if i == run else "cell removed" , s=1)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'R{name}_{date}_run{i}-{i+1}.png'))
        # plt.savefig(os.path.join(Plots, f'{date}', f'{name}_{date}_run{i}.png'))
        plt.show()

    def XY_vs_nu(self, lambda_path, lockins_path, name, run, n, B, power):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockins_path)

        if name == 'f':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, X1f, Y1f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_{\omega}$ ($\mu$V)', r'$\text{XY}_{\omega}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={-B} G, $P$={power} $\mu$W' + ' @'+ str(date))
        elif name == '2f':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, X2f, Y2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_{2\omega}$ ($\mu$V)', r'$\text{XY}_{2\omega}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={-B} G, $P$={power} $\mu$W' + ' @'+ str(date))
        elif name == 'dc':
            self.XYplot(Bristol_t, Lambda, para, lockins_t, Xdc, Ydc, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{XY}_\text{dc}$ (mV)', r'$\text{XY}_\text{dc}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={-B} G, $P$={power} $\mu$W' + ' @'+ str(date))
            
    def R_vs_nu(self, lambda_path, lockins_path, name, run, n, B, power):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.FR_double_Kvapor(lockins_path)

        if name == 'f':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, R1f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_{\omega}$ ($\mu$V)', r'$\text{R}_{\omega}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={-B} G, $P$={power} $\mu$W' + ' @'+ str(date))
        elif name == '2f':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, R2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_{2\omega}$ ($\mu$V)', r'$\text{R}_{2\omega}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={-B} G, $P$={power} $\mu$W' + ' @'+ str(date))
        elif name == 'dc':
            self.Rplot(Bristol_t, Lambda, para, lockins_t, Rdc, run-1, n, name, 'Frequency (GHz)',
                                  r'$\text{R}_\text{dc}$ (mV)', r'$\text{R}_\text{dc}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={-B} G, $P$={power} $\mu$W' + ' @'+ str(date))
            
plotter = Plot()
date_input = '05-09-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.XY_vs_nu(Bristol_path, Lockins_path, 'dc', 11, 5, 5.12, 4.99) 
# plotter.R_vs_nu(Bristol_path, Lockins_path, 'f', 5, 5, 5.07, 3) 