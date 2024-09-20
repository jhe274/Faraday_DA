import os, glob
import datetime as dt
import numpy as np
import scipy.special
from scipy.optimize import curve_fit
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
    
    def find_closest_indices(self, arrays, power, n):
        arrays = [np.asarray(subarray, dtype=np.float64) for subarray in arrays]
        diffs = [np.abs(subarray - power) for subarray in arrays]
        
        # Create a structured array to keep track of subarray index, element index, value, and difference
        all_closest = np.array([(sub_idx, idx, subarray[idx], diff[idx])
                                for sub_idx, (subarray, diff) in enumerate(zip(arrays, diffs))
                                for idx in range(len(subarray))],
                            dtype=[('sub_idx', int), ('idx', int), ('value', float), ('diff', float)])
        
        # Sort by the difference and select the top `n` closest values
        sorted_closest = np.sort(all_closest, order='diff')[:n]
        
        # Extract the result as a list of tuples
        result = [(item['sub_idx'], item['idx'], item['value']) for item in sorted_closest]
        return result
    
    def corr_plot(self, dates, phytype, datatype, power, n):
        
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        results = [self.reader.ellip_theta(glob.glob(os.path.join(processed_path, date, '*.csv'))) for date in dates]
        T, B, P, x1, y1, y2 = zip(*[(res[1], res[2], res[3], res[4], res[5], res[6]) for res in results])

        # Convert to float arrays
        T, B, P = map(self.analyzer.convert_to_float_array, [T, B, P])
        
        # Find closest indices
        result = self.find_closest_indices(P, power, n)

        freq, theta = [], []
        for sub_idx, idx, value in result:
            nu = self.consts.c / x1[sub_idx][idx]                                                                                       # [Hz]
            detun = (nu - self.consts.Nu39_D2) * 1e-6                                                                                   # [GHz]
            label = f'T={T[sub_idx][idx]:.2f}°C, $B_z$={B[sub_idx][idx]:.2f} G, P={P[sub_idx][idx]:.2f} μW'

            if phytype == 'CD':
                if datatype == 'raw':
                    ax.plot(detun, y1[sub_idx][idx] * 1e6, '.', label=label, markersize=4)
                # elif datatype == 'grad':
            elif phytype == 'CB':
                if datatype == 'raw':
                    ax.plot(detun, y2[sub_idx][idx] * 1e6, '.', label=label, markersize=4)
                elif datatype == 'grad':
                    freq.append(nu)
                    theta.append(y2[sub_idx][idx])
        
        Delta_the_6G, Delta_the_5G = [], []
        for idx, x_val in enumerate(freq[2]):
            detun = (freq[2] - self.consts.Nu39_D2) * 1e-6 
            idx_6G = np.argmin(np.abs(freq[1] - x_val))
            Delta_the_6G.append(theta[1][idx_6G] - theta[2][idx])
        grad_6G = np.gradient(Delta_the_6G, .98)
        ax.plot(detun, grad_6G * 1e3, '.', label='6G-5G', markersize=4)
        for idx, x_val in enumerate(freq[0]):
            detun = (freq[0] - self.consts.Nu39_D2) * 1e-6 
            idx_5G =  np.argmin(np.abs(freq[2] - x_val))
            Delta_the_5G.append(theta[2][idx_5G] - theta[0][idx])
        grad_5G = np.gradient(Delta_the_5G, 1.02)
        ax.plot(detun, grad_5G * 1e3, '.', label='5G-4G', markersize=4)

        self.plot_settings(power, phytype, datatype)
    
    def plot_settings(self, power, phytype, datatype):
        plt.xlabel(r'Frequency (MHz)', fontsize=25)
        plt.xticks(np.arange(-5000, 6000, 100), fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(-1000,1000)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        plt.legend(loc='best', markerscale=10, fontsize=25)

        plots_path = os.path.join(processed_path, 'Correlation plots', 'Sort by power', f'{power:.0f} Microwatts')
        os.makedirs(plots_path, exist_ok=True)

        if phytype == 'CD':
            if datatype == 'raw':
                plt.ylabel(r'Ellipticity (μrad.)', fontsize=25)
                plt.title(f'Ellipticity vs Frequency Detuning with $P\simeq${power:.0f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_Ellipticity_vs_Detuning_{power:.0f}microW.png'))
            elif datatype == 'grad':
                plt.ylabel(r'$\nabla_{\nu}\epsilon$ (μrad/mG)', fontsize=25)
                plt.title(rf'$\nabla_\nu\epsilon$ vs Frequency Detuning with $P\simeq${power:.0f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_gradEllipticity_vs_Detuning_{power:.0f}microW.png'))
        elif phytype == 'CB':
            if datatype == 'raw':
                plt.ylabel(r'Faraday Rotation (μrad.)', fontsize=25)
                plt.title(f'Faraday Rotation vs Frequency Detuning with $P\simeq${power:.2f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_FR_vs_Detuning_{power:.0f}microW(zoomed).png'))
            elif datatype == 'grad':
                plt.ylabel(r'$\nabla_B\theta$ (μrad/mG)', fontsize=25)
                plt.title(rf'$\nabla_B\theta$ vs Frequency Detuning with $P\simeq${power:.0f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_gradFR_vs_Detuning_{power:.0f}microW.png'))
        plt.show()

if __name__ == "__main__":
    dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    # dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed data')

    plotter = Plot()
    dates = ['06-18-2024', '06-24-2024', '06-25-2024', '06-26-2024', '06-27-2024', '07-01-2024',
                '05-07-2024', '05-09-2024', '05-15-2024', '05-19-2024',
                '05-23-2024', '05-29-2024', '05-31-2024', '06-05-2024', '06-07-2024']
    plotter.corr_plot(dates, 'CB', 'grad', 400, 3)