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

    def process_and_plot(self, Bristol_t, Lambda, para, y_t, R, run, n, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        for i in range(run, run+2):
        # i = run
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], R[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], y_t[i], R[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], y_t[i], para[i][2], n)
            Lambd, y = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R[i][l_idx])
            x = self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9                                                   # Frequency: [GHz]
            # y = y * 1e3                                                                                                     # RMS Voltage: [mV]                            
            y = y * 1e6                                                                                                     # RMS Voltage: [microV]

            ax.scatter(x, y[1:], label=f'{r"L->H" if i == run else r"H->L"} Wide Scan', s=5)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'{name}_{date}_run{i}-{i+1}.png'))
        # plt.savefig(os.path.join(Plots, f'{date}', f'{name}_{date}_run{i}.png'))
        plt.show()


    def y_vs_nu(self, lambda_path, lockins_path, name, run, n, B, power):
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.Double_modu_theta(lockins_path)

        if name == 'R1f':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, R1f, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{1f}$ ($\mu$V)', r'$R_{1f}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} nW' + ' @'+ str(date))
        elif name == 'R2f':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, R2f, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{2f}$ ($\mu$V)', r'$R_{2f}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} nW' + ' @'+ str(date))
        elif name == 'Rdc':
            self.process_and_plot(Bristol_t, Lambda, para, lockins_t, Rdc, run-1, n, name, 'Frequency (GHz)',
                                  r'$R_{dc}$ (mV)', r'$R_{dc}$ vs Frequency, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} nW' + ' @'+ str(date))

plotter = Plot()
date_input = '03-19-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.y_vs_nu(Bristol_path, Lockins_path, 'R1f', 7, 5, 5.103, 873) 