import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

fig, ax = plt.subplots()
cwd = os.getcwd()
dir_path = os.path.join(cwd, 'Research', 'PhD Project', 'Faraday Rotation Measurements', 'K vapor cell')
Bristol = os.path.join(dir_path, 'Bristol data')
Lockins = os.path.join(dir_path, 'Lock-ins data')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()

    def theta_mea_vs_nu(self, lambda_path, lock_in_path, i, n):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, theta = self.Analyze.Double_mod_theta(lock_in_path)
        # print(lockins_t[0][0], lockins_t[0][-1])
        # para, Lock_in_t, theta = Analyze.Triple_mod_theta(lock_in_path)
        Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
        Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i])
        l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], Lock_in_t[i], para[i][2], n)
        Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
        Thet = Thet[l_idx[1:-50]]

        print(np.where(np.isnan(Lambd[:-50])))
        # print(len(Lambd[:-50]), len(Thet[Thet != 0]))
        nu0 = self.Consts.c / self.Consts.Lambda_D2 * 1e-9 # [GHz]
        x = self.Consts.c / Lambd[:-50] * 1e-9 - nu0 # [MHz]
        y = Thet[Thet != 0] * 1e3 # [millirad]

        x1 = x[:np.argmin(x)]
        x2 = x[np.argmin(x):]
        y1 = y[:np.argmin(x)]
        y2 = y[np.argmin(x):]
        valleys, _ = find_peaks(-y, prominence=.5, height=(-45, -30))
        ax.plot(x1, y1, color='black', label=r'forward scan')
        ax.plot(x2, y2, color='blue', label=r'backward scan')

        for j, k in enumerate(valleys):
            color = plt.cm.tab10(j)
            ax.scatter(x[k], y[k], color=color, s=20)
            ax.vlines(x=x[k], ymin=5, ymax=y[k], 
                      linestyle='--', color=color, label=r'$\nu$={:.6f} GHz, $\theta$={:.2f} millirad.'.format(x[k], y[k]), linewidth=2)
        
        plt.xlabel(r'Frequency (GHz)')
        plt.ylabel(r'Polarization Rotation (millirad.)')
        plt.xticks(np.arange(self.Consts.c/766.705824 - nu0, self.Consts.c/766.696020 - nu0, .5))
        plt.yticks(np.arange(5,50,5))
        # plt.xlim(Consts.c/766.706, Consts.c/766.696)
        plt.ylim(5, 45)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Frequency, $B$=11.92 G, $P$=1.142 mW @{date}')
        plt.show()

plotter = Plot()
# date_input = input("Enter the date (MM-DD-YYYY): ")
date_input = '01-24-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
# TC300_path = glob.glob(os.path.join(TC300, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.theta_mea_vs_nu(Bristol_path, Lockins_path, 8, 5)