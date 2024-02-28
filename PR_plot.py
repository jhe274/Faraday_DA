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

    def frequency_detuning(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.Analyze.Double_modu_theta(lockin_path)
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))

        for i in range(run-1, run+1):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x = self.Consts.c / Lambd * 1e-9 - self.Consts.Nu39_D2 * 1e-9                                                   # [GHz]
            y = Thet[1:] * 1e3                                                                                              # [millirad]

            ax.plot(x, y, label=f'{r"L->H" if i == run-1 else r"H->L"} Wide Scan')
        
        lambda_theo = np.linspace(766.69*1e-9, 766.71*1e-9, 2000)
        y_theo = self.Theory.FR_theta(lambda_theo, 0.0718, B*1e-4, 26, self.Consts.Lambda39_D1, self.Consts.Lambda39_D2)
        x_theo = self.Consts.c / lambda_theo * 1e-9 - self.Consts.Nu39_D2  * 1e-9
        plt.plot(x_theo, y_theo * 1e3, '--', color='red', label='Theory')

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.ylabel(r'Polarization Rotation (millirad.)', fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(-4,6)
        plt.ylim(0,2)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(f'Polarization Rotation vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'PR_vs_Frequency(theory)_{date}_run{run}-{run+1}.png'))
        plt.show()

    def theory_plot(self, l, B, T):
        x = np.linspace(760*1e-9, 777*1e-9, 20000)
        y = self.Theory.FR_theta(x, l, B*1e-4, T, self.Consts.Lambda39_D1, self.Consts.Lambda39_D2)

        plt.plot(x*1e9, y * 1e3, '--', color='red', label='Theory')
        plt.ylim(0, 2e-5)
        plt.grid(True)
        plt.show()

plotter = Plot()
date_input = '02-27-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.frequency_detuning(Bristol_path, Lockins_path, 9, 5, 5.103, 470)
# plotter.theory_plot(0.0718, 5.103, 26)