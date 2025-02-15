import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory_calculations import Theory
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze
from scipy.signal import find_peaks

class LaserDrift:
    def __init__(self):
        """
        Initialize the Plot class with constants, theory, reader, and analyzer modules.
        """
        self.consts = Consts()  # Load physical constants
        self.theory = Theory()  # Load theoretical models
        self.reader = Read()    # Utilities for reading data
        self.analyzer = Analyze()  # Data analysis utilities

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run)
    
    def process(self, t, y1, y2, dtype, run):
        if dtype == 'wavelength':
            scale_factor = 1e9  # wavelenght: [nm]
            y = y1 * scale_factor
            x = t / 60  # Convert time to minutes
            y_fit, slope, intercept, y_mean, residuals, y_std = self.analyzer.drift_fit(x, y)
        elif dtype == 'frequency':
            scale_factor = 1e-9  # frequency: [GHz]
            y = y2 * scale_factor
            x = t / 60  # Convert time to minutes
            y_fit, slope, intercept, y_mean, residuals, y_std = self.analyzer.drift_fit(x, y)

        return x, y, y_fit, slope, intercept, y_mean, residuals, y_std
    
    def plot_process(self, timestamp, wavelength, dtype, run, name, ave_unit, std_unit, xlabel, ylabel, title, date):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        t, wl, nu, x, y, y_fit, slope, y_mean, y_std = [], [], [], [], [], [], [], [], []
        for i in self.number_of_runs(run):
            t, wl = self.analyzer.filter_data(timestamp[i], wavelength[i])
            nu = np.array([self.consts.c / wl[j] for j in range(len(wl))])  # [GHz]
            x, y, y_fit, slope, intercept, y_mean, residuals, y_std = self.process(t, wl, nu, dtype, run)

            lower_bound = slope*x + intercept - y_std
            upper_bound = slope*x + intercept + y_std

            if i == run-1:
                unitfactor = 1e3
                ax.plot(x, y, color='C3', alpha=0.4, linestyle='-', linewidth=0.5, \
                        label=fr'$\overline{{{name}}}$={round(y_mean,3):.3f} {ave_unit}')
                ax.plot(x, y_fit, '--', color='r', 
                    label=fr'$\nabla_t{name}$={round(slope * unitfactor,3):.3f} {std_unit}/min')
                ax.fill_between(x, lower_bound, upper_bound, facecolor='C3', alpha=0.4, label=f'$\\sigma$={round(y_std * unitfactor,3):.3f} {std_unit}')
            else:
                unitfactor = 1e3
                ax.plot(x, y, color='C0', alpha=0.4, linestyle='-', linewidth=0.5, \
                        label=fr'$\overline{name}$={round(y_mean,3):.3f} {ave_unit}')
                ax.plot(x, y_fit, '--', color='b',  
                    label=fr'$\nabla_t{name}$={round(slope * unitfactor,3):.3f} {std_unit}/min')
                ax.fill_between(x, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(y_std * unitfactor,3):.3f} {std_unit}')

        ax.ticklabel_format(useOffset=False, style='plain')
        plt.xlabel(xlabel, fontsize=25)
        plt.xticks(fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(False)
        ax.legend(loc='best', fontsize=25)
        # plt.title(title, fontsize=25)
        plt.savefig(os.path.join(plots, f'{date}', f'{dtype}_vs_time_{date}_run{run}.png'))
        plt.show()
        
    def wavelength_frequency(self, wavelengthmeter_path, date, run, dtype):
        # Read wavelength data from Bristol wavelength meter
        timestamp, wavelength = self.reader.read_bristol(wavelengthmeter_path)

        if dtype == 'wavelength':
            name = r'\lambda'
            ave_unit = r'nm'
            std_unit = r'pm'
            self.plot_process(timestamp, wavelength, dtype, run, name, ave_unit, std_unit, 
                            r'Time (min)', r'Wavelength (nm)', f'Wavelength vs Time, run{run}' + ' @'+ str(date), date)
        elif dtype == 'frequency':
            name = r'\nu'
            ave_unit = r'GHz'
            std_unit = r'MHz'
            self.plot_process(timestamp, wavelength, dtype, run, name, ave_unit, std_unit,  
                            r'Time (min)', r'Frequency (GHz)', f'Frequency vs Time, run{run}' + ' @'+ str(date), date)

if __name__ == "__main__":
    dir_path = os.path.join(
    os.path.expanduser('~'),  # Directory path on personal computer
    'OneDrive', 
    'Files', 
    'Graduate_study', 
    'Research', 
    'PhD_project', 
    'Faraday_rotation_measurements'
    )
    # dir_path = os.path.join(os.getcwd(),  # Directory path on Faraday lab computer
    # 'Faraday_rotation_measurements', 
    # )
    K_vapor = os.path.join(dir_path, 'K_vapor_cell')
    wavelengthmeter = os.path.join(K_vapor, 'Wavelengthmeter_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data')
    
    plotter = LaserDrift()
    date_input = '02-15-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    wavelengthmeter_path = glob.glob(os.path.join(wavelengthmeter, date, '*.csv'))
    plotter.wavelength_frequency(wavelengthmeter_path, date, 2, 'frequency')