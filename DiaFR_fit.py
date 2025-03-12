import os, glob
import datetime as dt
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory_calculations import Theory
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

class Plot:
    """
    A class to handle data visualization, analysis, and processing for Faraday rotation measurements.
    """
    def __init__(self):
        """
        Initialize the Plot class with constants, theory, reader, and analyzer modules.
        """
        self.consts = Consts()  # Load physical constants
        self.theory = Theory()  # Load theoretical models
        self.reader = Read()    # Utilities for reading data
        self.analyzer = Analyze()  # Data analysis utilities
        self.theory = Theory()
        self.l = (7.5 - 0.159 * 2) * 1e-2  # Optical path length in meters

    def convert_to_float(self, data_tuple):
        """Convert a tuple of lists of strings into a tuple of lists of floats."""
        return tuple([np.array([float(item.strip()) for item in sublist]) for sublist in data_tuple])
    
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
    
    def plot_number_density_fit(self, dates, power, n):
        # Read the data from the specified path
        results = [self.reader.read_processed_data(glob.glob(os.path.join(processed_path, date, '*.csv'))) for date in dates]
        temps, Bzs, powers, x1, y1, y2 = zip(*[(res[1], res[2], res[3], res[4], res[5], res[6]) for res in results])

        # Convert to float
        temps, Bzs, powers = map(self.convert_to_float, [temps, Bzs, powers])

        # Find closest indices
        result = self.find_closest_indices(powers, power, n)

        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        temp, Bz, power, freq, detune, theta = [], [], [], [], [], []
        for sub_idx, idx, value in result:
            nu = self.consts.c / x1[sub_idx][idx][::10] # [Hz]
            detuning = (nu - self.consts.K39_D2_Hz) * 1e-9 # [GHz]

            temp.append(temps[sub_idx][idx]) # [°C]
            Bz.append(Bzs[sub_idx][idx]) # [G]
            power.append(powers[sub_idx][idx]) # [µW]
            freq.append(nu) # [Hz]
            detune.append(detuning) # [GHz]
            theta.append(y2[sub_idx][idx][::10] * 2) # [rad]

        color = ['b', 'g', 'r']
        for i in range(len(temp)):
            label = f'T={temp[i]:.2f}°C, $B_z$={Bz[i]:.2f} G, power={power[i]:.2f} μW'
            ax.plot(detune[i], theta[i] * 1e6, '.', color=color[i], alpha=0.4, markersize=10, label=label)

            # Fit the data
            initial_guess = [[1, 20.3, -6.11, -0.002, 1000, 1000, 60], [1, 21.85, -4.05, -0.002, 1000, 1000, 60], [1, 20.7, -5.12, -0.002, 1000, 1000, 60]]
            bounds = [([1, 20., -6.115, -.1, 0, 0, -100], [5, 21, -6.105, .1, 2000, 4000, 100]), 
                      ([1, 21.5, -4.055, -.1, 0, 0, -100], [5, 22, -4.045, .1, 2000, 4000, 100]), 
                        ([1, 20.5, -5.125, -.1, 0, 0, -100], [3, 21, -5.115, .1, 2000, 4000, 100])]
            fitted_y, fit_Kn, fit_temp, fit_Bz, fit_PK, fit_gamma1, fit_gamma2, fit_offset = self.number_density_fit(freq[i], theta[i], initial_guess[i], bounds[i])
            ax.plot(detune[i], fitted_y * 1e6, '--', color=color[i], label=rf'[K]={fit_Kn:.2f}$\times10^{{8}}$ cm$^{{-3}}$, $B_z$={fit_Bz:.2f} G, P={fit_PK*1e2:.2f}%', linewidth=2)

        ax.set_xlabel(r'Frequency Detuning, $\nu$ (GHz)', fontsize=25)
        ax.set_ylabel(r'Faraday rotation, $\theta$ (µrad)', fontsize=25)
        ax.set_xticks(np.arange(-5, 6, 1))
        ax.tick_params(axis='both', which='major', labelsize=25)
        ax.legend(loc='best', fontsize=25)
        save_path = os.path.join(plots, 'Atomic_density_fit.png')
        plt.savefig(save_path)
        plt.show()
        
    def number_density_fit(self, x, y, initial_guess, bounds):
        params, covariance = curve_fit(self.theory.resonant_FR, x, y, p0=initial_guess, bounds=bounds)
        fitted_y = self.theory.resonant_FR(x, params[0], params[1], params[2], params[3], params[4], params[5], params[6])

        return fitted_y, params[0], params[1], params[2], params[3], params[4], params[5], params[6]


if __name__ == "__main__":
    dir_path = os.path.join(
    os.path.expanduser('~'),  # Expands to your home directory
    'OneDrive', 
    'Files', 
    'Graduate_study', 
    'Research', 
    'PhD_project', 
    'Faraday_rotation_measurements'
    )
    # dir_path = os.path.join(os.getcwd(),  
    # 'Faraday_rotation_measurements', 
    # )
    K_vapor = os.path.join(dir_path, 'K_vapor_cell')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data')
    
    plotter = Plot()
    dates = ['05-07-2024', '06-07-2024', '06-18-2024']
    plotter.plot_number_density_fit(dates, 0.5, 3)


