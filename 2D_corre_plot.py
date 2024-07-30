import os, glob
import datetime as dt
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()
    
    def corr_plot(self, path, phytype):
        
        date, temp, Bz, power, wl, ellip, theta = self.reader.ellip_theta(path)
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        x = []
        
        for i in range(len(date)):
            x.append(self.consts.c / wl[i] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                   # [GHz]
            if phytype == 'CD':
                ax.plot(x[i], ellip[i]*1e6, '.', label=f'T={temp[i]}°C, $B_z$={Bz[i]} G, P={power[i]} $\mu$W', markersize=2)
                
            if phytype == 'CB':
                ax.plot(x[i], theta[i]*1e6, '.', label=f'T={temp[i]}°C, $B_z$={Bz[i]} G, P={power[i]} $\mu$W', markersize=2)

        self.plot_settings(date[0], phytype)
        
    def plot_settings(self, date, phytype):
        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(-5,6)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        plt.legend(loc='best', fontsize=25)
        if phytype == 'CD':
            plt.ylabel(r'Ellipticity (microrad.)', fontsize=25)
            plt.title(f'Ellipticity  vs Frequency @{date}', fontsize=25)
            plt.savefig(os.path.join(Plots, f'[X]_Ellipticity_vs_Frequency_{date}.png'))
        elif phytype == 'CB':
            plt.ylabel(r'Faraday Rotation (microrad.)', fontsize=25)
            plt.title(f'Faraday Rotation vs Frequency @{date}', fontsize=25)
            plt.savefig(os.path.join(Plots, f'[X]_FR_vs_Frequency_{date}.png'))
        plt.show()

if __name__ == "__main__":
    dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    # dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed data')

    plotter = Plot()
    dates_4G = ['06-18-2024', '06-24-2024', '06-25-2024', '06-26-2024', '06-27-2024', '07-01-2024']
    dates_5G = ['05-07-2024', '05-09-2024', '05-15-2024', '05-19-2024', '05-23-2024']
    dates_6G = ['05-29-2024', '05-31-2024', '06-05-2024', '06-07-2024']
    date_input = '05-19-2024'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    data_path = glob.glob(os.path.join(processed_path, date, '*.csv'))
    Plots = os.path.join(processed_path, date, 'Plots')
    plotter.corr_plot(data_path, 'CB')