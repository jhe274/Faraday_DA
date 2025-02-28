from scipy.fft import fft, fftfreq
import numpy as np
import matplotlib.pyplot as plt

import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory_calculations import Theory
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run)

    def plot_process(self, x, y, component, run, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in self.number_of_runs(run):

            if name == '1f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == '2f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'dc':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'm2f':
                scale_factor = 1e3  # RMS Voltage: [mV]

            if component == 'X':
                label = (r'$\text{X}_\text{f}$' if name == '1f' else 
                        r'$\text{X}_\text{2f}$' if name == '2f' else 
                        r'$\text{X}_\text{dc}$' if name == 'dc' else
                        r'$\text{X}_\text{m2f}$')
            elif component == 'Y':
                label = (r'$\text{Y}_\text{f}$' if name == '1f' else
                        r'$\text{Y}_\text{2f}$' if name == '2f' else
                        r'$\text{Y}_\text{dc}$' if name == 'dc' else
                        r'$\text{Y}_\text{m2f}$')
            elif component == 'R':
                label = (r'$\text{R}_\text{1f}$' if name == '1f' else
                        r'$\text{R}_\text{2f}$' if name == '2f' else
                        r'$\text{R}_\text{dc}$' if name == 'dc' else
                        r'$\text{R}_\text{m2f}$')
            if i == run-1:
                # Apply the scaling factor to X and Y arrays
                y[i] = y[i] * scale_factor

                # Recompute the FFT
                N = len(y[i])  # Number of data points
                sampling_interval = np.median(np.diff(x[i]))  # Average time step in seconds

                # Compute FFT and corresponding frequencies
                fft_values = fft(y[i] - np.mean(y[i]))  # Remove DC component
                frequencies_fft = fftfreq(N, d=sampling_interval)  # Frequency axis

                # Convert to positive frequencies only
                positive_freqs = frequencies_fft[:N // 2]
                fft_magnitudes = np.abs(fft_values[:N // 2])

                # Convert frequency to period (in minutes) while avoiding division by zero
                valid_indices = positive_freqs > 0
                positive_freqs = positive_freqs[valid_indices]
                fft_magnitudes = fft_magnitudes[valid_indices]
                periods_fft_minutes = 1 / positive_freqs / 60  # Convert to minutes

                ax.plot(periods_fft_minutes, fft_magnitudes, label=label)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        # plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(1, max(periods_fft_minutes)) 
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(plots, f'{date}', rf'FFT_of_{component}{name}, run' + f'{run}.png'))
        plt.show()

    def XYR_vs_time(self, lockins_path, name, component, run):
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f = self.reader.read_lockins(lockins_path)
        para, lockins_t, R1f, R2f, Rdc, Rm2f = self.analyzer.R_lockins(lockins_path)

        if component == 'X':
            if name == '1f':
                self.plot_process(lockins_t, X1f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of X1f, run' + f'{run}')
            elif name == '2f':
                self.plot_process(lockins_t, X2f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of X2f, run' + f'{run}')
            elif name == 'dc':
                self.plot_process(lockins_t, Xdc, component, run, name, '"Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Xdc, run' + f'{run}')
            elif name == 'm2f':
                self.plot_process(lockins_t, Xm2f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Xm2f, run' + f'{run}')
        elif component == 'Y':
            if name == '1f':
                self.plot_process(lockins_t, Y1f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Y1f, run' + f'{run}')
            elif name == '2f':
                self.plot_process(lockins_t, Y2f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Y2f, run' + f'{run}')
            elif name == 'dc':
                self.plot_process(lockins_t, Ydc, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Ydc, run' + f'{run}')
            elif name == 'm2f':
                self.plot_process(lockins_t, Ym2f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Ym2f, run' + f'{run}')
        elif component == 'R':
            if name == '1f':
                self.plot_process(lockins_t, R1f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of R1f, run' + f'{run}')
            elif name == '2f':
                self.plot_process(lockins_t, R2f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of R2f, run' + f'{run}')
            elif name == 'dc':
                self.plot_process(lockins_t, Rdc, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Rdc, run' + f'{run}')
            elif name == 'm2f':
                self.plot_process(lockins_t, Rm2f, component, run, name, 'Period (minutes)',
                                    r'FFT Magnitude', rf'FFT spectrum of Rm2f, run' + f'{run}')
    
if __name__ == "__main__":
    dir_path = os.path.join(
    # os.path.expanduser('~'),  # Directory path on personal computer
    'D:',
    'OneDrive', 
    'Files', 
    'Graduate_study', 
    'Research', 
    'PhD_project', 
    'Faraday_rotation_measurements'
    )
    print(dir_path)
    # dir_path = os.path.join(os.getcwd(),  # Directory path on Faraday lab computer
    # 'Faraday_rotation_measurements', 
    # )
    K_vapor = os.path.join(dir_path, 'K_vapor_cell')
    lockins = os.path.join(K_vapor, 'Lockins_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '02-27-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    lockins_path = glob.glob(os.path.join(lockins, date, '*.lvm'))
    # for i in range(1,8):
    plotter.XYR_vs_time(lockins_path, '2f', 'X', 6)