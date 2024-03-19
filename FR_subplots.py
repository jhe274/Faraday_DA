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
    
    def y_vs_Frequency(self, lambda_path, lockin_path, run, n):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.Double_modu_theta(lockin_path)
        x, y = [], []
        fig, ax = plt.subplots(2, 1, figsize=(25, 12))  # 4 subplots, each with 1 column

        for i in range(run-1, run+3):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            # Bristol_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], epsilon[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            # Lambd, Epsi = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            Lambd, Thet = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                             # [GHz]
            # y.append(Epsi[1:] * 1e3)                                                                                        # [millirad]
            y.append(Thet[1:] * 1e3)                                                                                        # [millirad]

        ax[0].plot(x[0], y[0], label=f'{"L->H"} Wide Scan')
        ax[0].plot(x[1], y[1], label=f'{"H->L"} Wide Scan')
        ax[0].tick_params(axis='both', labelsize=15)
        ax[0].set_xticks(np.arange(-5, 6, 1))
        ax[0].set_title(r'Vapor cell removed, $P$=860 nW, $\lambda_{PEM}$=766.7 nm, $A$=2.391 rad', fontsize=20)
        ax[0].grid(True)
        ax[0].legend(loc='best', fontsize=15)

        ax[1].plot(x[2], y[2], label=f'{"L->H"} Wide Scan')
        ax[1].plot(x[3], y[3], label=f'{"H->L"} Wide Scan')
        ax[1].tick_params(axis='both', labelsize=15)
        ax[1].set_xticks(np.arange(-5, 6, 1))
        ax[1].set_title(r'Vapor cell inserted, $P$=860 nW, $\lambda_{PEM}$=766.7 nm, $A$=2.391 rad', fontsize=20)
        ax[1].grid(True)
        ax[1].legend(loc='best', fontsize=15)

        # ax[1, 0].plot(x[4], y[4], label=f'{"L->H"} Wide Scan')
        # ax[1, 0].plot(x[5], y[5], label=f'{"H->L"} Wide Scan')
        # ax[1, 0].tick_params(axis='both', labelsize=15)
        # ax[1, 0].set_title(r'Vapor cell removed, $P$=550 nW, $\lambda_{PEM}$=766.691 nm', fontsize=20)
        # ax[1, 0].grid(True)
        # ax[1, 0].legend(loc='best', fontsize=15)

        # ax[1, 1].plot(x[6], y[6], label=f'{"L->H"} Wide Scan')
        # ax[1, 1].plot(x[7], y[7], label=f'{"H->L"} Wide Scan')
        # ax[1, 1].tick_params(axis='both', labelsize=15)
        # # ax[1, 1].set_xlabel('Frequency (GHz)', fontsize=20)
        # # ax[1, 1].set_ylabel('Polarization Rotation (millirad.)', fontsize=20)
        # ax[1, 1].set_title(r'Vapor cell inserted, $P$=570 nW, $\lambda_{PEM}$=766.691 nm', fontsize=20)
        # ax[1, 1].grid(True)
        # ax[1, 1].legend(loc='best', fontsize=15)
            
        fig.text(0.5, 0.04, 'Frequency (GHz)', ha='center', va='center', fontsize=25)
        # fig.text(0.02, 0.5, 'Ellipticity (millirad.)', ha='center', va='center', rotation='vertical', fontsize=25)
        fig.text(0.02, 0.5, 'Faraday Rotation (millirad.)', ha='center', va='center', rotation='vertical', fontsize=25)
        # plt.suptitle('Ellipticity measurements at $B_z$=5.103 G, $T=$25°C', fontsize=35)
        plt.suptitle('Faraday rotation measurements at $B_z$=5.103 G, $T=$25°C', fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        # plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity_subplots_@{date}.png'))
        plt.savefig(os.path.join(Plots, f'{date}', f'FR_subplots_@{date}.png'))
        plt.show()

plotter = Plot()
date_input = '03-18-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.y_vs_Frequency(Bristol_path, Lockins_path, 1, 5)