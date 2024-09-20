import os, glob
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze
from scipy.signal import find_peaks

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()
    
    def corr_plot(self, dates, phytype, datatype):
        
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        
        for i in dates:
            path = glob.glob(os.path.join(processed_path, i, '*.csv'))
            date, T, B, P, x1, y1, y2 = self.reader.ellip_theta(path)
            T, B, P = map(self.analyzer.convert_to_float_array, [T, B, P])
            x = []
            for j in range(len(path)):
                x.append(self.consts.c / x1[j] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                   # [GHz]
                smoothed_y1 = self.analyzer.smooth(y1[j]*1e6, 90)
                smoothed_y2 = self.analyzer.smooth(y2[j]*1e6, 90)
                y1_grad = np.gradient(y1[j]*1e6, x[j])
                y2_grad = np.gradient(y2[j]*1e6, x[j])
                label = f'P={P[j]:.2f} μW'
                if phytype == 'CD':
                    if datatype == 'raw':
                        ax.plot(x[j], y1[j]*1e6, '.', label=label, markersize=2)
                        ax.plot(x[j], smoothed_y1, '.', label=label, markersize=2)
                    if datatype == 'grad':
                        ax.plot(x[j], y1_grad, '.', label=label, markersize=2)
                elif phytype == 'CB':
                    if datatype == 'raw':
                        ax.plot(x[j], y2[j]*1e6, '.', label=label, markersize=2)
                        ax.text(-3, 1700, r'$B_z\simeq4.0$ G', fontsize=70, color='black', ha='center', va='bottom')
                    if datatype == 'grad':
                        ax.plot(x[j], y2_grad, '.', label=label, markersize=2)

        self.plot_settings(phytype, datatype)

    def plot_settings(self, phytype, datatype):
        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(-5,6)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        plt.legend(loc='best', markerscale=10, fontsize=15)

        plots_path = os.path.join(processed_path, 'Correlation plots', 'Sort by B-field', '6 Gauss measurements')
        os.makedirs(plots_path, exist_ok=True)

        if phytype == 'CD':
            if datatype == 'raw':
                plt.ylabel(r'Ellipticity (microrad.)', fontsize=25)
                plt.title(r'Ellipticity vs Frequency Detuning with $B_z\simeq6.1$ G', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_Ellipticity_vs_Detuning_6G.png'))
            elif datatype == 'grad':
                plt.ylabel(r'$\nabla_\nu\epsilon$ (microrad/GHz)', fontsize=25)
                plt.title(r'$\nabla_\nu\epsilon$ vs Frequency Detuning with $B_z\simeq6.1$ G', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_gradEllipticity_vs_Detuning_6G.png'))
        elif phytype == 'CB':
            if datatype == 'raw':
                plt.ylabel(r'Faraday Rotation (microrad.)', fontsize=25)
                plt.title(r'Faraday Rotation vs Frequency Detuning with $B_z\simeq4.0$ G', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_FR_vs_Detuning_4G(poster).png'))
            elif datatype == 'grad':
                plt.ylabel(r'$\nabla_\nu\theta$ (microrad/GHz)', fontsize=25)
                plt.title(r'$\nabla_\nu\theta$ vs Frequency Detuning with $B_z\simeq6.1$ G', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_gradFR_vs_Detuning_6G.png'))
        plt.show()

if __name__ == "__main__":
    dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    # dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed data')

    plotter = Plot()
    dates_4G = ['06-18-2024', '06-24-2024', '06-25-2024', '06-26-2024', '06-27-2024', '07-01-2024']
    dates_5G = ['05-07-2024', '05-09-2024', '05-15-2024', '05-19-2024']
    dates_6G = ['05-23-2024', '05-29-2024', '05-31-2024', '06-05-2024', '06-07-2024']
    
    plotter.corr_plot(dates_4G, 'CB', 'raw')