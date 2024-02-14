import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

fig, ax = plt.subplots()
dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
# dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
K_vapor = os.path.join(dir_path, 'K vapor cell')
Bristol = os.path.join(K_vapor, 'Bristol data')
Lockins = os.path.join(K_vapor, 'Lockins data')
DLCpro = os.path.join(K_vapor, 'TopticaDLCpro data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()
    
    def frequency_detuning(self, lambda_path, lockin_path, run, n):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, R2f, Rdc, theta = self.Analyze.Double_modu_theta(lockin_path)
        x, y = [], []
        fig, ax = plt.subplots(1, 2, figsize=(35, 12))  # 2 subplots, each with 1 column

        for i in range(run-1, run+11):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x.append(self.Consts.c / Lambd * 1e-9 - self.Consts.Nu39_D2 * 1e-9)                             # [GHz]
            y.append(Thet[1:] * 1e3)                                                                      # [millirad]

        # ax[0, 0].plot(x[0], y[0], label=f'{"L->H"} Wide Scan')
        # ax[0, 0].plot(x[1], y[1], label=f'{"H->L"} Wide Scan')
        # ax[0, 0].tick_params(axis='x', labelsize=15)
        # ax[0, 0].tick_params(axis='y', labelsize=15)
        # ax[0, 0].set_xlabel('Frequency (GHz)', fontsize=20)
        # ax[0, 0].set_ylabel('Polarization Rotation (millirad.)', fontsize=20)
        # ax[0, 0].set_title('$P=18.78 \mu$W', fontsize=20)
        # ax[0, 0].grid(True)
        # ax[0, 0].legend(loc='best')

        # ax[0, 1].plot(x[2], y[2], label=f'{"L->H"} Wide Scan')
        # ax[0, 1].plot(x[3], y[3], label=f'{"H->L"} Wide Scan')
        # ax[0, 1].tick_params(axis='x', labelsize=15)
        # ax[0, 1].tick_params(axis='y', labelsize=15)
        # ax[0, 1].set_xlabel('Frequency (GHz)', fontsize=20)
        # ax[0, 1].set_ylabel('Polarization Rotation (millirad.)', fontsize=20)
        # ax[0, 1].set_title('$P=15.04 \mu$W', fontsize=20)
        # ax[0, 1].grid(True)
        # ax[0, 1].legend(loc='best')

        # ax[1, 0].plot(x[4], y[4], label=f'{"L->H"} Wide Scan')
        # ax[1, 0].plot(x[5], y[5], label=f'{"H->L"} Wide Scan')
        # ax[1, 0].tick_params(axis='x', labelsize=15)
        # ax[1, 0].tick_params(axis='y', labelsize=15)
        # ax[1, 0].set_xlabel('Frequency (GHz)', fontsize=20)
        # ax[1, 0].set_ylabel('Polarization Rotation (millirad.)', fontsize=20)
        # ax[1, 0].set_title('$P=10.00 \mu$W', fontsize=20)
        # ax[1, 0].grid(True)
        # ax[1, 0].legend(loc='best')

        # ax[1, 1].plot(x[6], y[6], label=f'{"L->H"} Wide Scan')
        # ax[1, 1].plot(x[7], y[7], label=f'{"H->L"} Wide Scan')
        # ax[1, 1].tick_params(axis='x', labelsize=15)
        # ax[1, 1].tick_params(axis='y', labelsize=15)
        # ax[1, 1].set_xlabel('Frequency (GHz)', fontsize=20)
        # ax[1, 1].set_ylabel('Polarization Rotation (millirad.)', fontsize=20)
        # ax[1, 1].set_title('$P=7.35 \mu$W', fontsize=20)
        # ax[1, 1].grid(True)
        # ax[1, 1].legend(loc='best')

        ax[0].plot(x[8], y[8], label=f'{"L->H"} Wide Scan')
        ax[0].plot(x[9], y[9], label=f'{"H->L"} Wide Scan')
        ax[0].tick_params(axis='x', labelsize=25)
        ax[0].tick_params(axis='y', labelsize=25)
        ax[0].set_title('$P=4.48 \mu$W', fontsize=25)
        ax[0].grid(True)
        ax[0].legend(loc='best')

        ax[1].plot(x[10], y[10], label=f'{"L->H"} Wide Scan')
        ax[1].plot(x[11], y[11], label=f'{"H->L"} Wide Scan')
        ax[1].tick_params(axis='x', labelsize=25)
        ax[1].tick_params(axis='y', labelsize=25)
        ax[1].set_title('$P=2.28 \mu$W', fontsize=25)
        ax[1].grid(True)
        ax[1].legend(loc='best')
            
        fig.text(0.5, 0.04, 'Detuning (GHz)', ha='center', va='center', fontsize=25)
        fig.text(0.02, 0.5, 'Polarization Rotation (millirad.)', ha='center', va='center', rotation='vertical', fontsize=25)
        plt.suptitle('Polarization rotation measurements at $B=$5.17 G, $T=$23°C', fontsize=35)
        # plt.tight_layout()
        plt.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
        plt.savefig(os.path.join(Plots, f'{date}', f'PR_subplots_@{date}.png'))

plotter = Plot()
date_input = '01-31-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.frequency_detuning(Bristol_path, Lockins_path, 1, 5)