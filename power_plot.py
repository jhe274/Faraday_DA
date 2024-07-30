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
    
    def convert_to_float_array(self, arr):
        return np.array([np.array(lst, dtype=np.float64) for lst in np.array(arr, dtype=object)], dtype=object)
    
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
    
    def corr_plot(self, dates, phytype, power, n):
        
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        results = [self.reader.ellip_theta(glob.glob(os.path.join(processed_path, date, '*.csv'))) for date in dates]
        T, B, P, x1, y1, y2 = zip(*[(res[1], res[2], res[3], res[4], res[5], res[6]) for res in results])

        # Convert to float arrays
        T, B, P = map(self.convert_to_float_array, [T, B, P])
        
        # Find closest indices
        result = self.find_closest_indices(P, power, n)

        for sub_idx, idx, value in result:
            detun = self.consts.c / x1[sub_idx][idx] * 1e-9 - self.consts.Nu39_D2 * 1e-9
            label = f'T={T[sub_idx][idx]:.2f}°C, $B_z$={B[sub_idx][idx]:.2f} G, P={P[sub_idx][idx]:.2f} $\mu$W'
            
            if phytype == 'CD':
                ax.plot(detun, y1[sub_idx][idx] * 1e6, '.', label=label, markersize=4)
            elif phytype == 'CB':
                ax.plot(detun, y2[sub_idx][idx] * 1e6, '.', label=label, markersize=4)
        
        self.plot_settings(power, phytype)

    def plot_settings(self, power, phytype):
        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(-5,6)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        plt.legend(loc='best', markerscale=5, fontsize=15)

        plots_path = os.path.join(processed_path, 'Correlation plots', 'Sort by power', f'{power:.0f} Microwatts')
        os.makedirs(plots_path, exist_ok=True)

        if phytype == 'CD':
            plt.ylabel(r'Ellipticity (microrad.)', fontsize=25)
            plt.title(f'Ellipticity vs Frequency Detuning with $P\simeq${power:.0f} $\mu$W', fontsize=25)
            plt.savefig(os.path.join(plots_path, f'[X]_Ellipticity_vs_Detuning_{power:.0f}microW.png'))
        elif phytype == 'CB':
            plt.ylabel(r'Faraday Rotation (microrad.)', fontsize=25)
            plt.title(f'Faraday Rotation vs Frequency Detuning with $P\simeq${power:.0f} $\mu$W', fontsize=25)
            plt.savefig(os.path.join(plots_path, f'[X]_FR_vs_Detuning_{power:.0f}microW.png'))
        plt.show()

if __name__ == "__main__":
    dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    # dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed data')

    plotter = Plot()
    dates = ['06-18-2024', '06-24-2024', '06-25-2024', '06-26-2024', '06-27-2024', '07-01-2024',
                '05-07-2024', '05-09-2024', '05-15-2024', '05-19-2024',
                '05-23-2024', '05-29-2024', '05-31-2024', '06-05-2024', '06-07-2024']
    plotter.corr_plot(dates, 'CB', 500, 3)