import os, glob
import datetime as dt
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
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
        self.l = (7.5 - 0.159 * 2) * 1e-2  # Optical path length in meters

    def convert_to_float(self, data_tuple):
        """Convert a tuple of lists of strings into a tuple of lists of floats."""
        return tuple([np.array([float(item.strip()) for item in sublist]) for sublist in data_tuple])
    
    def find_closest_indices(self, arrays, parameter, n):
        arrays = [np.asarray(subarray, dtype=np.float64) for subarray in arrays]
        diffs = [np.abs(subarray - parameter) for subarray in arrays]
        
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
            step = 1
            nu = self.consts.c / x1[sub_idx][idx][::step] # [Hz]
            detuning = (nu - self.consts.K39_D2_Hz) * 1e-9 # [GHz]

            temp.append(temps[sub_idx][idx]) # [°C]
            Bz.append(Bzs[sub_idx][idx]) # [G]
            power.append(powers[sub_idx][idx]) # [µW]
            freq.append(nu) # [Hz]
            detune.append(detuning) # [GHz]
            theta.append(y2[sub_idx][idx][::step] * 2) # [rad]

        color = ['r', 'g', 'b']
        lines = []
        for i in [1, 0, 2]:
            label = f'$B_z$={Bz[i]:.2f} G, $T$={temp[i]:.2f} °C'
            line, = ax.plot(detune[i], theta[i] * 1e3, '.', c=color[i], alpha=1, markersize=2, label=label)
            lines.append(line)  # Store handles for legend

            # # Fit the data
            # initial_guess = [[1, 21.85, -4.05, -0.002, 1000, 1000, 60], [1, 20.7, -5.12, -0.002, 1000, 1000, 60], [1, 20.3, -6.11, -0.002, 1000, 1000, 60]]
            # bounds = [([1, 21.5, -4.055, -.1, 0, 0, -100], [5, 22, -4.045, .1, 2000, 4000, 100]), 
            #             ([1, 20.5, -5.125, -.1, 0, 0, -100], [3, 21, -5.115, .1, 2000, 4000, 100]), 
            #             ([1, 20., -6.115, -.1, 0, 0, -100], [5, 21, -6.105, .1, 2000, 4000, 100]), ]
            # fitted_y, fit_Kn, fit_temp, fit_Bz, fit_PK, fit_gamma1, fit_gamma2, fit_offset = self.number_density_fit(freq[i], theta[i], initial_guess[i], bounds[i])
            # ax.plot(detune[i], fitted_y * 1e6, '--', color=color[i], label=rf'[K]={fit_Kn:.2f}$\times10^{{8}}$ cm$^{{-3}}$', linewidth=2)

        ax.set_xlabel(r'Frequency Detuning, $\nu$ (GHz)', fontsize=30)
        ax.set_ylabel(r'Faraday Rotation, $\theta$ (rad)', fontsize=30)
        ax.set_xticks(np.arange(-5, 6, 1))
        ax.tick_params(axis='both', which='major', labelsize=30)
        # ax.legend(loc='best', fontsize=30)

        # Define a custom legend handler that increases only the marker size
        class HandlerLineMarker(HandlerLine2D):
            def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                # Increase the marker size for the legend without changing the plot markers
                line = super().create_artists(legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans)
                for l in line:
                    l.set_markersize(10)  # Set a larger marker size in the legend
                return line

        # Customizing legend marker size
        ax.legend(handles=lines, fontsize=30, loc="best",
                handler_map={line: HandlerLineMarker() for line in lines})  # Use custom handler

        plt.tight_layout()
        save_path = os.path.join(plots, 'Atomic_density_fit(tightlayout).png')
        plt.savefig(save_path)
        plt.show()
        
    def number_density_fit(self, x, y, initial_guess, bounds):
        params, covariance = curve_fit(self.theory.resonant_FR, x, y, p0=initial_guess, bounds=bounds)
        fitted_y = self.theory.resonant_FR(x, params[0], params[1], params[2], params[3], params[4], params[5], params[6])

        return fitted_y, params[0], params[1], params[2], params[3], params[4], params[5], params[6]
    
    def plot_power_dependence(self, dates, Bz, n):
        # Read the data from the specified path
        results = [self.reader.read_processed_data(glob.glob(os.path.join(processed_path, date, '*.csv'))) for date in dates]
        temps, Bzs, powers, x1, y1, y2 = zip(*[(res[1], res[2], res[3], res[4], res[5], res[6]) for res in results])

        # Convert to float
        temps, Bzs, powers = map(self.convert_to_float, [temps, Bzs, powers])

        # Find closest indices
        result = self.find_closest_indices(Bzs, Bz, n)

        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        temp, Bz, power, freq, detune, theta = [], [], [], [], [], []
        for sub_idx, idx, value in result:
            step = 1
            nu = self.consts.c / x1[sub_idx][idx][::step]
            detuning = (nu - self.consts.K39_D2_Hz) * 1e-9

            temp.append(temps[sub_idx][idx]) # [°C]
            Bz.append(Bzs[sub_idx][idx]) # [G]
            power.append(powers[sub_idx][idx]) # [µW]
            freq.append(nu) # [Hz]
            detune.append(detuning) # [GHz]
            theta.append(y2[sub_idx][idx][::step] * 2) # [rad]

        # indices_4G = [6, 9, 10, 14, 16, 18, 19, 20, 22, 24]
        indices_6G = [6, 9, 10, 14, 16, 18, 19, 20, 22, 23]
        lines = []
        for i in indices_6G:
            intensity = power[i] * 1e-3 / (np.pi * (0.22/2) ** 2)
            label = f'$I$={intensity:.3f} mW/cm$^2$'
            label = f'$P$={power[i]:.1f} µW'
            line, = ax.plot(detune[i], theta[i] * 1e3, '.', alpha=1, markersize=2, label=label)
            lines.append(line)  # Store handles for legend

        ax.set_xlabel(r'Frequency Detuning, $\nu$ (GHz)', fontsize=30)
        ax.set_ylabel(r'Faraday Rotation, $\theta$ (rad)', fontsize=30)
        ax.set_xticks(np.arange(-5, 6, 1))
        ax.tick_params(axis='both', which='major', labelsize=30)
        # ax.legend(loc='best', fontsize=30)

        # Define a custom legend handler that increases only the marker size
        class HandlerLineMarker(HandlerLine2D):
            def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                # Increase the marker size for the legend without changing the plot markers
                line = super().create_artists(legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans)
                for l in line:
                    l.set_markersize(20)  # Set a larger marker size in the legend
                return line

        # Customizing legend marker size
        ax.legend(handles=lines, fontsize=30, loc="best",
                handler_map={line: HandlerLineMarker() for line in lines})  # Use custom handler

        plt.tight_layout()
        save_path = os.path.join(plots, 'FR_power_dependence_6G(tightlayout).png')
        plt.savefig(save_path)
        plt.show()
        

    def plot_slope(self, dates, power, n):
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
            step = 1
            nu = self.consts.c / x1[sub_idx][idx][::step] # [Hz]
            detuning = (nu - self.consts.K39_D2_Hz) * 1e-9 # [GHz]

            temp.append(temps[sub_idx][idx]) # [°C]
            Bz.append(Bzs[sub_idx][idx]) # [G]
            power.append(powers[sub_idx][idx]) # [µW]
            freq.append(nu) # [Hz]
            detune.append(detuning) # [GHz]
            theta.append(y2[sub_idx][idx][::step] * 2) # [rad]

        color = ['r', 'g', 'b']
        lines = []
        for i in [1, 0]:
            label = f'$B_z$={Bz[i]:.2f} G, $T$={temp[i]:.2f} °C'
            line, = ax.plot(detune[i], theta[i] * 1e3 / np.abs(Bz[i]), '.', c=color[i], alpha=1, markersize=2, label=label)
            lines.append(line)  # Store handles for legend

        ax.set_xlabel(r'Frequency Detuning, $\nu$ (GHz)', fontsize=30)
        ax.set_ylabel(r'Slope ($10^{-3}$ rad/G)', fontsize=30)
        ax.set_xticks(np.arange(-5, 6, 1))
        ax.tick_params(axis='both', which='major', labelsize=30)
        # ax.legend(loc='best', fontsize=30)

        # Define a custom legend handler that increases only the marker size
        class HandlerLineMarker(HandlerLine2D):
            def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                # Increase the marker size for the legend without changing the plot markers
                line = super().create_artists(legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans)
                for l in line:
                    l.set_markersize(10)  # Set a larger marker size in the legend
                return line

        # Customizing legend marker size
        ax.legend(handles=lines, fontsize=30, loc="best",
                handler_map={line: HandlerLineMarker() for line in lines})  # Use custom handler

        plt.tight_layout()
        save_path = os.path.join(plots, 'Slope_300microW_2inputs(tightlayout).png')
        plt.savefig(save_path)
        plt.show()

    def plot_sensitivity(self, dates, power, n):
        # Read the data from the specified path
        results = [self.reader.read_processed_data(glob.glob(os.path.join(processed_path, date, '*.csv'))) for date in dates]
        temps, Bzs, powers, x1, y1, y2 = zip(*[(res[1], res[2], res[3], res[4], res[5], res[6]) for res in results])

        # Convert to float
        temps, Bzs, powers = map(self.convert_to_float, [temps, Bzs, powers])

        # Find closest indices
        result = self.find_closest_indices(powers, power, n)

        # Shot-noise-limited angular sensitivity
        spectral_responsivity = 0.56 # A/W
        QE = spectral_responsivity * self.consts.h * self.consts.c / (self.consts.e * 766.7e-9) # Quantum efficiency
        snl_sensitivity = np.sqrt(self.consts.h * self.consts.c / (2 * QE * 766.7e-9 * 54.52e-6)) # [rad/sqrt(Hz)]
        print('Shot-noise-limited angular sensitivity = {:.2f} x 1e-8 G/sqrt(Hz)'.format(snl_sensitivity*1e8))

        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        temp, Bz, power, freq, detune, theta = [], [], [], [], [], []
        for sub_idx, idx, value in result:
            step = 1
            nu = self.consts.c / x1[sub_idx][idx][::step] # [Hz]
            detuning = (nu - self.consts.K39_D2_Hz) * 1e-9 # [GHz]

            temp.append(temps[sub_idx][idx]) # [°C]
            Bz.append(Bzs[sub_idx][idx]) # [G]
            power.append(powers[sub_idx][idx]) # [µW]
            freq.append(nu) # [Hz]
            detune.append(detuning) # [GHz]
            theta.append(y2[sub_idx][idx][::step] * 2) # [rad]

        color = ['r', 'C0', 'b']
        lines = []
        for i in [1]: # Plot with 4 G measurement
            label = f'$B_z$={Bz[i]:.2f} G, $T$={temp[i]:.2f} °C'
            line, = ax.plot(detune[i], np.abs(Bz[i]) * snl_sensitivity * 1e4 / theta[i], '.', c=color[i], alpha=1, markersize=5, label=label)
            lines.append(line)  # Store handles for legend
            print("Minimum ∂B_z/∂θ = {:.2f} G/rad".format(min(np.abs(Bz[i]/theta[i]))))
        ax.set_xlabel(r'Frequency Detuning, $\nu$ (GHz)', fontsize=30)
        ax.set_ylabel(r'Sensitivity, $\delta B_z$ ($10^{-4}$ G/$\sqrt{\text{Hz}}$)', fontsize=30)
        # ax.set_xticks(np.arange(-5, 6, 1))
        ax.set_xlim(-1.5,1.5)
        ax.set_ylim(0, 20)
        ax.tick_params(axis='both', which='major', labelsize=30)
        # ax.legend(loc='best', fontsize=30)

        # Define a custom legend handler that increases only the marker size
        class HandlerLineMarker(HandlerLine2D):
            def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                # Increase the marker size for the legend without changing the plot markers
                line = super().create_artists(legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans)
                for l in line:
                    l.set_markersize(10)  # Set a larger marker size in the legend
                return line

        # Customizing legend marker size
        # ax.legend(handles=lines, fontsize=30, loc="best",
                # handler_map={line: HandlerLineMarker() for line in lines})  # Use custom handler

        plt.tight_layout()
        save_path = os.path.join(plots, 'Sensitivity_300microW_2inputs(tightlayout).png')
        # plt.savefig(save_path)
        plt.show()

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
    dates = ['06-18-2024', '06-24-2024', '06-25-2024', '06-26-2024', '06-27-2024', '07-01-2024',
                '05-07-2024', '05-09-2024', '05-15-2024', '05-19-2024',
                '05-23-2024', '05-29-2024', '05-31-2024', '06-05-2024', '06-07-2024']
    # plotter.plot_number_density_fit(dates, 0.5, 3)
    plotter.plot_power_dependence(dates, 6, 24)
    # plotter.plot_slope(dates, 300, 3)
    # plotter.plot_sensitivity(dates, 300, 3)


