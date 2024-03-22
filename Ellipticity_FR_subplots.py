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
DLCpro = os.path.join(K_vapor, 'TopticaDLCpro data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()
    
    def y_vs_Frequency(self, lambda_path, lockin_path, run, n, B, P):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.FR_double_Kvapor(lockin_path)
        x, y = [], []
        fig, ax = plt.subplots(2, 2, figsize=(25, 12))  # 4 subplots, each with 1 column

        for i in range(run-1, run+3):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], epsilon[i])
            # Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Epsi = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            # Lambd, Thet = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                             # [GHz]
            y.append(Epsi[1:] * 1e3)                                                                                        # [millirad]
            # y.append(Thet[1:] * 1e3)                                                                                        # [millirad]

        ax[0, 0].plot(x[0], y[0], label=r'$\theta_p$=89.5300°, $\theta_a$=46.2435°')
        # ax[0, 0].plot(x[1], y[1], label=f'{"H->L"} Wide Scan')
        ax[0, 0].tick_params(axis='both', labelsize=15)
        ax[0, 0].set_xticks(np.arange(-5, 6, 1))
        ax[0, 0].set_title(f'Vapor cell removed, run{run}', fontsize=20)
        ax[0, 0].grid(True)
        ax[0, 0].legend(loc='best', fontsize=15)

        ax[0, 1].plot(x[1], y[1], label=r'$\theta_p$=89.5070°, $\theta_a$=46.2365°')
        # ax[0, 1].plot(x[3], y[3], label=f'{"H->L"} Wide Scan')
        ax[0, 1].tick_params(axis='both', labelsize=15)
        ax[0, 1].set_xticks(np.arange(-5, 6, 1))
        ax[0, 1].set_title(f'Vapor cell removed, run{run+1}', fontsize=20)
        ax[0, 1].grid(True)
        ax[0, 1].legend(loc='best', fontsize=15)

        ax[1, 0].plot(x[2], y[2], label=r'$\theta_p$=89.4835°, $\theta_a$=46.2290°')
        # ax[1, 0].plot(x[5], y[5], label=f'{"H->L"} Wide Scan')
        ax[1, 0].tick_params(axis='both', labelsize=15)
        ax[1, 0].set_xticks(np.arange(-5, 6, 1))
        ax[1, 0].set_title(f'Vapor cell removed, run{run+2}', fontsize=20)
        ax[1, 0].grid(True)
        ax[1, 0].legend(loc='best', fontsize=15)

        ax[1, 1].plot(x[3], y[3], label=r'$\theta_p$=89.4630°, $\theta_a$=46.2435°')
        # ax[1, 1].plot(x[7], y[7], label=f'{"H->L"} Wide Scan')
        ax[1, 1].tick_params(axis='both', labelsize=15)
        ax[1, 1].set_xticks(np.arange(-5, 6, 1))
        # ax[1, 1].set_xlabel('Frequency (GHz)', fontsize=20)
        # ax[1, 1].set_ylabel('Polarization Rotation (millirad.)', fontsize=20)
        ax[1, 1].set_title(f'Vapor cell removed, run{run+3}', fontsize=20)
        ax[1, 1].grid(True)
        ax[1, 1].legend(loc='best', fontsize=15)
            
        fig.text(0.5, 0.04, 'Frequency (GHz)', ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, 'Ellipticity (millirad.)', ha='center', va='center', rotation='vertical', fontsize=25)
        # fig.text(0.02, 0.5, 'Faraday Rotation (millirad.)', ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle(r'Ellipticity measurements at $B_z$=' + f'{B} G, ' + r'$P$=' + f'{P} nW, ' + r'$T=$25°C', fontsize=35)
        # plt.suptitle(r'Faraday rotation measurements at $B_z$=' + f'{B} G, ' + r'$P$=' + f'{P} nW, ' + r'$T=$25°C', fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity_subplots_@{date}.png'))
        # plt.savefig(os.path.join(Plots, f'{date}', f'FR_subplots_@{date}.png'))
        plt.show()

plotter = Plot()
date_input = '03-21-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.y_vs_Frequency(Bristol_path, Lockins_path, 1, 5, 5.103, 865)