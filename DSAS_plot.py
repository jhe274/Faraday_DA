import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze
import numpy.polynomial.polynomial as poly

fig, ax = plt.subplots()
# dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
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

    def process_and_plot(self, Bristol_t, Lambda, y_t, data, run, n):
        x, y = [], []
        for i in range(run, run+2):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], y_t[i], data[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], y_t[i], data[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], y_t[i], 0, n)
            x.append(self.Consts.c / Lambda[i][b_idx] * 1e-9 - self.Consts.Nu39_D2 * 1e-9)              # [GHz]
            y.append(data[i][l_idx])
        
        data = np.array([x[1], y[1]]).T
        data = data[np.argsort(data[:, 0])]
        X, Y = data[:, 0], data[:, 1]

        # Find indices of duplicate elements in X
        unique_X, unique_indices = np.unique(X, return_index=True)
        duplicate_indices = np.setdiff1d(np.arange(len(X)), unique_indices)

        # Remove duplicate elements from X
        X = np.delete(X, duplicate_indices)

        # Remove corresponding elements from Y
        Y = np.delete(Y, duplicate_indices)

        Y_cs = CubicSpline(X, Y)
        XX = np.linspace(X.min(), X.max(), 1000)
        YY = Y_cs(XX)

        savgol_params = {'window_length': 10, 'polyorder': 3}
        dY = np.gradient(YY, XX)
        dYs = savgol_filter(dY, deriv=1, delta=XX[1]-XX[0], **savgol_params)

        peaks, _ = find_peaks(Y, prominence=.005, height=(0, .03))
        valleys, _ = find_peaks(-Y, prominence=0.02, height=(0, 0.04))
        print(peaks, valleys)
        K39_diff = (X[valleys[0]] + self.Consts.Nu39_D2 - self.Consts.Nu39_D2_B) * 1e-9
        K41_D2_A = (self.Consts.Nu41_D2_A-self.Consts.Nu39_D2) * 1e-9 + K39_diff
        K41_D2_B = (self.Consts.Nu41_D2_B-self.Consts.Nu39_D2) * 1e-9 + K39_diff
        K41_D2_C = (self.Consts.Nu41_D2_C-self.Consts.Nu39_D2) * 1e-9 + K39_diff
        

        fig, axs = plt.subplots(2, 1, figsize=(12, 12))
        axs[0].plot(X, Y, color='red')
        for j,k in enumerate(peaks):
            axs[0].axvline(x=X[k], linestyle='--', color='red')
            axs[0].axvline(x=X[k], linestyle='--', color='red')
        for j,k in enumerate(valleys):
            axs[0].axvline(x=X[k], linestyle='--', color='red')
        axs[0].axvline(x=K41_D2_A, color='orange', linestyle='--')
        axs[0].axvline(x=K41_D2_B, color='orange', linestyle='--')
        axs[0].axvline(x=K41_D2_C, color='orange', linestyle='--')
        axs[0].grid(True)
        axs[0].set_xlim(-.8, 1.2)
        # axs[0].xaxis.set_major_formatter(ScalarFormatter(useMathText=True, useOffset=False))
        axs[0].tick_params(axis='x', labelsize=15)
        axs[0].tick_params(axis='y', labelsize=15)
        axs[0].set_ylabel('Intensity (au)', fontsize=20)
        axs[0].set_xlabel('Frequency (GHz)', fontsize=20)
        
        axs[1].plot(XX, dY, color='blue')
        # axs[1].plot(XX, dYs, color='blue')
        for j,k in enumerate(peaks):
            axs[1].axvline(x=X[k], linestyle='--', color='red')
            axs[1].axvline(x=X[k], linestyle='--', color='red')
        for j,k in enumerate(valleys):
            axs[1].axvline(x=X[k], linestyle='--', color='red')
        axs[1].grid(True)
        axs[1].set_xlim(-.8, 1.2)
        # axs[1].set_ylim(-2.2,1.5)
        axs[1].tick_params(axis='x', labelsize=15)
        axs[1].tick_params(axis='y', labelsize=15)
        axs[1].set_ylabel('Intensity (au)', fontsize=20)
        axs[1].set_xlabel('Frequency (GHz)', fontsize=20)
        plt.savefig(os.path.join(Plots, f'{date}', f'SAS_KVC_{date}_run{run+2}.png'))

    def V_vs_nu(self, lambda_path, y_path, run, n):
        """
        Plot function of Lock-in RMS values vs Wavelength
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        x, y, Y, DLCpro_t = self.Read.DLCpro_WideScan(y_path)
        self.process_and_plot(Bristol_t, Lambda, DLCpro_t, y, run-1, n)
                            # r'Doppler-free SAS of K vapor cell @$T=36°$C')
                            # r'Doppler-free SAS of K vapor cell with blanced phoodetector @$T=36°$C')
plotter = Plot()
# date_input = input("Enter the date (MM-DD-YYYY): ")
date_input = '02-12-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
DLCpro_path = glob.glob(os.path.join(DLCpro, date, '*.csv'))
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.V_vs_nu(Bristol_path, DLCpro_path, 3, 5)