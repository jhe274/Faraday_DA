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

    def Ellipticity_vs_Frequency(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.Double_modu_theta(lockin_path)
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))

        for i in range(run-1, run+1):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], epsilon[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Epsi = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            x = self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9                                                   # [GHz]
            y = Epsi[1:] * 1e3                                                                                              # [millirad]

            ax.plot(x, y, label=f'{r"L->H" if i == run-1 else r"H->L"} Wide Scan')

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.ylabel(r'Ellipticity (millirad.)', fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # plt.xlim(-4,6)
        # plt.ylim(0,2)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(f'Ellipticity vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity_vs_Frequency_{date}_run{run}-{run+1}.png'))
        plt.show()

    def FR_vs_Frequency(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.Double_modu_theta(lockin_path)
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))

        for i in range(run-1, run+1):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x = self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9                                                   # [GHz]
            y = Thet[1:] * 1e3                                                                                              # [millirad]

            ax.plot(x, y, label=f'{r"L->H" if i == run-1 else r"H->L"} Wide Scan')
        
        # lambda_theo = np.linspace(766.69*1e-9, 766.71*1e-9, 2000)
        # y_theo = self.Theory.FR_theta1(lambda_theo, 0.0718, B*1e-4, 26, self.Consts.Lambda39_D1, self.Consts.Lambda39_D2)
        # x_theo = self.Consts.c / lambda_theo * 1e-9 - self.Consts.Nu39_D2  * 1e-9
        # plt.plot(x_theo, y_theo * 1e3, '--', color='red', label='Theory')

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.ylabel(r'Faraday Rotation (millirad.)', fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # plt.xlim(-4,6)
        # plt.ylim(0,2)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(f'Faraday Rotation vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'PR_vs_Frequency_{date}_run{run}-{run+1}.png'))
        plt.show()

    def theory_plot(self, l, B, T):
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        x1 = np.linspace(760*1e-9, 777*1e-9, 20000)
        y1 = self.theory.FR_theta1(x1, l, B*1e-4, T, self.consts.Lambda39_D1, self.consts.Lambda39_D2)
        # plt.plot(x1 * 1e9, y1 * 1e3, '-', color='red', label='Theory1')

        x2 = np.linspace(self.consts.c/(777*1e-9), self.consts.c/(760*1e-9), 20000)
        y2 = self.theory.FR_theta2(x2, l, B*1e-4, T, self.consts.Nu39_D1, self.consts.Nu39_D2)
        plt.plot(x2 * 1e-9, y2 * 1e3, '--', color='blue', label='Theory1')
        plt.xticks(np.arange(min(x2*1e-9), max(x2*1e-9), (max(x2*1e-9)-min(x2*1e-9))/10))
        ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        # plt.ylim(0, 2e5)
        plt.grid(True)
        plt.show()

plotter = Plot()
date_input = '03-19-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
# plotter.Ellipticity_vs_Frequency(Bristol_path, Lockins_path, 7, 5, 5.103, 873)
plotter.FR_vs_Frequency(Bristol_path, Lockins_path, 7, 5, 5.103, 873)
# plotter.theory_plot(0.0718, 5.103, 26)