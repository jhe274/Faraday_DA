import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import find_peaks
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

fig, ax = plt.subplots()
# dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements', 'K vapor cell')
dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements', 'K vapor cell')
Bristol = os.path.join(dir_path, 'Bristol data')
Lockins = os.path.join(dir_path, 'Lock-ins data')

class Plot:

    def __init__(self):
        self.Consts = Consts()
        self.Theory = Theory()
        self.Read = Read()
        self.Analyze = Analyze()

    def theta_mea_vs_nu(self, lambda_path, lock_in_path, n):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.Read.Bristol(lambda_path)
        para, lockins_t, theta = self.Analyze.Double_mod_theta(lock_in_path)
        xmin, xmax = float('-inf'), float('inf')
        for i in range(6,8):
            Bristol_t[i], Lambda[i] = self.Analyze.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.Analyze.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.Analyze.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.Analyze.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            nu0 = self.Consts.c / self.Consts.Lambda_D2 * 1e-9                              # [GHz]
            x = self.Consts.c / Lambd * 1e-9 - nu0                                          # [GHz]
            y = Thet[1:] * 1e3                                                              # [millirad]

            valleys, _ = find_peaks(-y, prominence=.5, height=(-28, -23))
            color = 'black' if i == 6 else 'blue'
            ax.plot(x, y, color=color, label=f'{r"L->H" if i == 6 else r"H->L"} Wide Scan')

            xmin = max(xmin, min(x))
            xmax = min(xmax, max(x))
        lambda_theo = np.linspace(766.695*1e-9, 766.705*1e-9, 2000)
        theta_theo = self.Theory.FR_theta(lambda_theo, 0.0718, 5*1e-4, 35, self.Consts.Lambda_D1, self.Consts.Lambda_D2)
        x_theo = self.Consts.c / lambda_theo * 1e-9 - nu0
        print(x_theo)
        ax.plot(x_theo, theta_theo * 1e3, '--', color='red', label='Theory')

        for j, k in enumerate(valleys):
            color = plt.cm.tab10(j)
            ax.scatter(x[k], y[k], color=color, s=20)
            ax.vlines(x=x[k], ymin=16, ymax=y[k], 
                      linestyle='--', color=color, label=r'$\nu$={:.6f} GHz, $\theta$={:.2f} millirad.'.format(x[k], y[k]), linewidth=2)
        
        plt.xlabel(r'Frequency (GHz)')
        plt.ylabel(r'Polarization Rotation (millirad.)')
        plt.xticks(np.arange(np.floor(xmin), np.ceil(xmax), 0.5))
        # plt.yticks(np.arange(5,50,5))
        plt.ylim(0, 22)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best')
        plt.title(f'Polarization Rotation vs Frequency, $B$=5.15 G, $P$=2.64 $\mu$W @{date}')
        plt.show()

plotter = Plot()
# date_input = input("Enter the date (MM-DD-YYYY): ")
date_input = '01-25-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.theta_mea_vs_nu(Bristol_path, Lockins_path, 5)