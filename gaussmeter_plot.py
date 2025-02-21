import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, butter, filtfilt
from lmfit import Model
from sklearn.linear_model import RANSACRegressor

def sine_wave(t, offset, amplitude, phase, freq):
    frequency = freq * 60
    return offset + amplitude * np.sin(2 * np.pi * frequency * t + phase)

class Plot:

    def __init__(self):
        self.reader = Read()
        self.analyzer = Analyze()

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run)
    
    def process_temperature(self, temps, run):
        T_mean = np.mean(temps[run])
        T_std = np.std(temps[run], ddof=1) / np.sqrt(len(temps[run]))

        return T_mean, T_std
    
    def plot_fft(self, x, y, run, xlabel, ylabel, file_name):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in self.number_of_runs(run):
            if i == run-1:
                freqs, fft_values, dominant_freq, amplitude_spectrum = self.analyzer.fft_peak(x[i], y[i])
                
                ax.plot(freqs, np.abs(fft_values)/len(y[i]), label='FFT', color='b')
                ax.set_xlim(0, 5)
                ax.set_xlabel(xlabel, fontsize=25)
                ax.set_ylabel(ylabel, fontsize=25)
                ax.tick_params(axis='x', labelsize=25)
                ax.tick_params(axis='y', labelsize=25)
        # plt.grid(False)
        # save_path = os.path.join(plots, date, file_name)
        # plt.savefig(save_path)
        plt.show()
    
    def noise_floor_plot(self, x, y, run, xlabel, ylabel, file_name):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in self.number_of_runs(run):
            if i == run-1:
                freqs, fft_values, dominant_freq, amplitude_spectrum = self.analyzer.fft_peak(x[i], y[i])
                noise_floor, psd, peak_power= self.analyzer.noise_floor(x[i], y[i])

                ax.plot(freqs[:len(y[i])//2], 10 * np.log10(psd[:len(y[i])//2]), color='C0', label='Power Spectral Density')
                ax.axhline(10 * np.log10(noise_floor), color='r', linestyle="--", label="Noise Floor")
                ax.scatter(dominant_freq, 10 * np.log10(peak_power), color='red', label="Signal Peak", zorder=3)
                ax.set_xlabel(xlabel, fontsize=25)
                ax.set_ylabel(ylabel, fontsize=25)
                ax.tick_params(axis='x', labelsize=25)
                ax.tick_params(axis='y', labelsize=25)
        plt.grid(False)
        plt.legend(loc='best', fontsize=25)
        save_path = os.path.join(plots, date, file_name)
        plt.savefig(save_path)
        plt.show()
    
    def bandpass_filter(self, data, lowcut, highcut, fs, order=4):
        nyquist = 0.5 * fs
        low = lowcut / nyquist
        high = highcut / nyquist
        b, a = butter(order, [low, high], btype='band', analog=False)
        
        return filtfilt(b, a, data)
        
    def sine_fit(self, x, y, y_mean, y_std, freq):
        """Curve fit using Scipy.optimize.curve_fit"""
        # initial guess
        initial_guess = [y_mean, 2*y_std, np.pi/2, freq]
        # Fit a sine wave to the data
        popt, pcov = curve_fit(sine_wave, x, y, p0=initial_guess)
        # popt, _ = curve_fit(self.sine_wave, x[i]/60, y1_filtered, p0=initial_guess)
        # Extract the fitted parameters
        fitted_offset, fitted_amplitude, fitted_phase, fitted_frequency = popt
        amplitude_uncertainty = np.sqrt(pcov[1, 1])
        # Print the true and fitted parameters
        print(f"Fitted Parameters:")
        print(f"Frequency = {fitted_frequency}")
        print(f"Offset = {fitted_offset}")
        print(f"Amplitude = {fitted_amplitude}")
        print(f"Phase = {fitted_phase}")
        # Generate the fitted sine wave
        y_sinefit = sine_wave(x, *popt)
        residuals = y - y_sinefit
        sigma = np.std(residuals, ddof=1)
        chi2_value = np.sum((residuals / sigma) ** 2)
        dof = len(y) - len(popt)
        print(f"Reduced Chi-squared: {chi2_value/dof:.3f}")

        return y_sinefit, fitted_offset, fitted_amplitude, fitted_phase, fitted_frequency, amplitude_uncertainty, chi2_value, dof

    def plot_process(self, x, y1, y2, run, xlabel, ylabel_y1, ylabel_y2, file_name):
        # Create a figure and axes for plotting
        fig, ax1 = plt.subplots(1, 1, figsize=(25.60, 14.40))
        # Create second Y-axis
        ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis
        
        for i in self.number_of_runs(run):
            if i == run-1:
                y1[i] = -y1[i]  # Invert the magnetic field data
                y1_linearfit, slope, intercept, y1_mean, residuals, y1_std = self.analyzer.drift_fit(x[i]/60, y1[i])
                T_mean, T_std = self.process_temperature(y2, i)
                freqs, fft_values, dominant_freq, amplitude_spectrum = self.analyzer.fft_peak(x[i], y1[i])

                # Apply band-pass filter around detected frequency
                lowcut, highcut = dominant_freq - 0.1, dominant_freq + 0.1
                y1_filtered = self.bandpass_filter(y1[i], lowcut, highcut, fs=10)

                y_sinefit, fitted_offset, fitted_amplitude, fitted_phase, fitted_frequency, amplitude_uncertainty, chi2_value, dof = self.sine_fit(x[i]/60, y1[i], y1_mean, y1_std, dominant_freq)

                lower_bound = fitted_offset - y1_std
                upper_bound = fitted_offset + y1_std

                ax1.plot(x[i]/60, y1[i], color='C0', alpha=0.4, label=f'$\\overline{{B_z}}$={round(y1_mean,3):.3f} G')
                # plot the linear fit
                # ax1.plot(x[i]/60, y1_linearfit, '--', color='b', label=fr'$\nabla_t B_z$={round(slope*1e3,1):.1f} mG/min')

                # plot the sine fit
                ax1.plot(x[i]/60, y_sinefit, '--', color='b', label=fr'$\Delta B_z$={round(2*fitted_amplitude*1e3,2):.2f}±{round(amplitude_uncertainty*1e3,2):.2f} mG')

                ax1.fill_between(x[i]/60, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(y1_std*1e3)} mG')
                # ax1.fill_between(x[i]/60, upper_bound, y1[i], where=y1[i] > upper_bound, fc='red', alpha=0.4, interpolate=True)
                # ax1.fill_between(x[i]/60, lower_bound, y1[i], where=y1[i] < lower_bound, fc='red', alpha=0.4, interpolate=True)

                ax1.set_xlabel(xlabel, fontsize=25)
                ax1.set_ylabel(ylabel_y1, color='C0', fontsize=25)
                ax1.tick_params(axis='x', labelsize=25)
                ax1.tick_params(axis='y', labelsize=25)

                ax2.plot(x[i]/60, y2[i], label=f'$\\overline{{T}}$={round(T_mean,2):.2f}±{round(T_std,2):.2f}°C', color='black', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=2000)
                ax2.set_ylabel(ylabel_y2, fontsize=25)
                # ax2.set_yticks(np.arange(22.06, 22.07, 0.002))
                ax2.tick_params(axis='y', labelsize=25)
            else:
                y1_linearfit, slope, intercept, y1_mean, residuals, y1_std = self.analyzer.drift_fit(x[i]/60, y1[i])
                T_mean, T_std = self.process_temperature(y2, i)

                lower_bound = slope*x[i]/60 + intercept - y1_std
                upper_bound = slope*x[i]/60 + intercept + y1_std

                ax1.plot(x[i]/60, y1[i], color='C3', alpha=0.4, label=f'$\\overline{{B_z}}$={round(y1_mean,3):.3f} G')
                ax1.plot(x[i]/60, y1_linearfit, '--', color='r', label=fr'$\nabla_t B_z$={round(slope*1e3,3):.3f} mG/min')
                ax1.fill_between(x[i]/60, lower_bound, upper_bound, facecolor='C3', alpha=0.4, label=f'$\\sigma$={round(y1_std*1e3)} mG')
                # ax1.fill_between(x[i]/60, upper_bound, y1[i], where=y1[i] > upper_bound, fc='r', alpha=0.4, interpolate=True)
                # ax1.fill_between(x[i]/60, lower_bound, y1[i], where=y1[i] < lower_bound, fc='r', alpha=0.4, interpolate=True)

                ax2.plot(x[i]/60, y2[i], label=f'$\\overline{{T}}$={round(T_mean,2):.2f}±{round(T_std,2):.2f}°C', color='black', linestyle='-', linewidth=1, 
                        marker='x', markersize=10, markevery=500)

            plt.title(f'$\chi^2/\\text{{dof}}$={round(chi2_value/dof,3):.3f}', fontsize=25)
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=20)
        plt.grid(False)
        save_path = os.path.join(plots, date, file_name)
        plt.savefig(save_path)
        plt.show()
    
    def gaussmter_vs_time(self, gaussmeter_path, run):
        timestamps, B0s, temps = self.reader.read_gaussmeter(gaussmeter_path)
        
        self.plot_process(timestamps, B0s, temps, run, 'Time (min)',
                            r'Magnetic flux density (G)', r'Temperature (°C)', 
                            f'Magnetic_field_and_temperature_{date}_run{run}.png')
        
        # self.plot_fft(timestamps, B0s, run, 'Frequency (Hz)', 'Amplitude', f'FFT_{date}_run{run}.png')

        # self.noise_floor_plot(timestamps, B0s, run, 'Frequency (Hz)', 'Power (dB)', f'Noise_floor_{date}_run{run}.png')
    

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
    gaussmeter = os.path.join(K_vapor, 'Gaussmeter_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '02-18-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    gaussmeter_path = glob.glob(os.path.join(gaussmeter, date, '*.csv'))
    plotter.gaussmter_vs_time(gaussmeter_path, 7)