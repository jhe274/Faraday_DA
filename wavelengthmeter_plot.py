import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from constants import Constants as Consts
from theory_calculations import Theory
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze
import allantools
from scipy.ndimage import uniform_filter1d

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
            x = t / 3600  # Convert time to minutes
            y_fit, slope, intercept, y_weightedmean, residuals, y_residualstd = self.analyzer.drift_fit(x, y, 100)
        elif dtype == 'frequency':
            scale_factor = 1e-9  # frequency: [GHz]
            y = y2 * scale_factor
            x = t / 60  # Convert time to minutes
            y_fit, slope, intercept, y_weightedmean, residuals, y_residualstd = self.analyzer.drift_fit(x, y, 100)

        return x, y, y_fit, slope, intercept, y_weightedmean, residuals, y_residualstd
    
    def plot_fft(self, timestamp, wavelength, run, xlabel, ylabel, file_name):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in self.number_of_runs(run):
            t, wl = self.analyzer.filter_data(timestamp[i], wavelength[i])
            wl_fit, slope, intercept, wl_weightedmean, residuals, wl_residualstd = self.analyzer.drift_fit(t, wl, 100)

            freqs, fft_values, amplitude_spectrum = self.analyzer.fft_peak(t, residuals)

            # Filter out negative frequencies
            positive_freqs = freqs[freqs > 0]
            positive_fft_values = fft_values[freqs > 0]

            if i == run - 1:
                ax.plot(positive_freqs, np.abs(positive_fft_values) / len(wl), color='b', label=f'run{run}')
            else:
                ax.plot(positive_freqs, np.abs(positive_fft_values) / len(wl), color='r', label=f'run{run+1}')

            ax.set_xlabel(xlabel, fontsize=25)
            ax.set_ylabel(ylabel, fontsize=25)
            ax.tick_params(axis='x', labelsize=25)
            ax.tick_params(axis='y', labelsize=25)

        # plt.grid(False)
        plt.legend(loc='best', fontsize=25)
        save_path = os.path.join(plots, date, file_name)
        plt.savefig(save_path)
        plt.show()
    
    def plot_process(self, timestamp, wavelength, dtype, run, name, ave_unit, std_unit, xlabel, ylabel, date):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        t, wl, nu, x, y, y_fit, slope, y_mean, y_std = [], [], [], [], [], [], [], [], []
        for i in self.number_of_runs(run):
            # Filter the data and convert wavelength to frequency
            t, wl = self.analyzer.filter_data(timestamp[i], wavelength[i])
            nu = np.array([self.consts.c / wl[j] for j in range(len(wl))])  # [GHz]
            instrument_uncertainty = 0.2 * 1e-6 * np.mean(wl[0]) # [m]
            sigma_f = self.consts.c * instrument_uncertainty * 1e-9 / (np.mean(wl[0])**2)  # [GHz]
            
            if i == run-1:
                unitfactor = 1e3 # frequency: [MHz]
                
                x1, y1, y1_fit, y1_slope, y1_intercept, y1_weightedmean, y1_residuals, y1_residualstd = self.process(t, wl, nu, dtype, run)
                y1_smoothed = uniform_filter1d(y1, size=int(0.32 * 100))  # τ = 0.3s, 100 Hz sampling
                # Fitting Scan with WideScan wavelength/frequency measurements
                # k0 = y1_slope # slope guess
                # b0 = min(y1) # intercept guess
                # a0 = 16 * 1e-3 # amplitude guess 7.8
                # phi0 = 0 # phase guess  
                # y1_scanfit, fitted_k1, fitted_b1, fitted_a1, fitted_phi1, sigma_a1, sigma_k1, sigma_a1, sigma_phi1, chi2_value1, dof1 = self.analyzer.drift_sine_fit(x1*60, y1, k0, b0, a0, phi0, sigma_f)
                
                # Fitting locked laser wavelength/frequency measurements
                k0 = y1_slope # slope guess
                b0 = min(y1) # intercept guess
                y1_linearfit, fitted_k1, fitted_b1, sigma_k1, sigma_b1, residuals, std, chi2_value1, dof1 = self.analyzer.linear_fit(x1, y1, y1_slope, y1_intercept, sigma_f)

                # Calculate the error bounds
                lower_bound1 = y1_linearfit - std
                upper_bound1 = y1_linearfit + std
                
                # Plot the measured data and fit for WideScan
                # ax.plot(x1, y1, color='C3', alpha=0.6, linestyle='-', linewidth=0.5, 
                #         label=fr'$\dot{name}$={round(fitted_k1*unitfactor,1):.1f} ± {round(sigma_drift*unitfactor,1):.1f} {std_unit}/s')
                # ax.plot(x1, y1_scanfit, '--', color='r', label=fr'$\Delta\nu$={round(2*np.abs(fitted_a1)*unitfactor,1):.1f} ± {round(sigma_a1*unitfactor,1):.1f} MHz')
                
                # Plot the measured data for locked laser
                # ax.plot(x1, y1, color='C3', alpha=0.6, linestyle='-', linewidth=0.5, \
                #         label=fr'$\overline{{{name}}}$={round(y1_weightedmean,3):.3f} {ave_unit}')
                
                # Plot the smoothened data
                ax.plot(x1, y1_smoothed, color='C3', alpha=1, linestyle='-', linewidth=1, \
                        label=fr'$\overline{{{name}}}$={round(y1_weightedmean,3):.3f} {ave_unit}')
                
                # Plot the linear fit
                ax.plot(x1, y1_linearfit, '--', color='r', 
                    label=fr'$\dot{name}$={round(fitted_k1 * unitfactor * 60,1):.1f} MHz/h')

                # Plot 1 sigma error band
                ax.fill_between(x1, lower_bound1, upper_bound1, facecolor='C3', alpha=0.4, label=f'$\\sigma$={round(std*unitfactor,1):.1f} {std_unit}')
                
                ax.ticklabel_format(useOffset=False, style='plain')
                ax.set_xlabel(xlabel, fontsize=25)
                plt.xticks(fontsize=25)
                ax.set_ylabel(ylabel, fontsize=25)
                plt.yticks(fontsize=25)
            else:
                unitfactor = 1e3 # frequency: [MHz]
                x2, y2, y2_fit, y2_slope, y2_intercept, y2_weightedmean, y2_residuals, y2_residualstd = self.process(t, wl, nu, dtype, run)
                
                # Fitting Scan with WideScan wavelength/frequency measurements
                k0 = y2_slope # slope guess
                b0 = min(y2) # intercept guess
                a0 = 16 * 1e-3 # amplitude guess
                phi0 = 0 # phase guess  
                y2_scanfit, fitted_slope2, fitted_intercept2, fitted_amplitude2, fitted_phase2, amplitude_uncertainty2, slope_uncertainty2, chi2_value2, dof2 = self.analyzer.drift_sine_fit(x2*60, y2, k0, b0, a0, phi0, sigma_f)
                
                # Calculate the drift uncertainty
                sem = y2_residualstd / np.sqrt(len(y2))
                sigma_drift= sigma_f / max(x2*60) # [GHz/s]
                
                # Calculate the error bounds
                lower_bound2 = y2_scanfit - y2_residualstd
                upper_bound2 = y2_scanfit + y2_residualstd

                # Plot measured data and fit for WideScan
                ax.plot(x2, y2, color='C0', alpha=0.6, linestyle='-', linewidth=0.5, 
                        label=fr'$\dot{name}$={round(fitted_slope2*unitfactor,2):.1f} ± {round(sigma_drift*unitfactor,1):.1f} {std_unit}/s')
                ax.plot(x2, y2_scanfit, '--', color='b', label=fr'$\Delta\nu$={round(2*np.abs(fitted_amplitude2)*unitfactor,1):.1f} ± {round(amplitude_uncertainty2*unitfactor,1):.1f} MHz')

                # Plot 1 sigma error band
                ax.fill_between(x2, lower_bound2, upper_bound2, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(y2_residualstd*unitfactor,1):.1f} {std_unit}')

                # --------- Add Inset Zoomed Plot ---------
                # axins = inset_axes(ax, width="50%", height="50%", 
                #                     bbox_to_anchor=(0.48, -0.2, 0.5, 0.5),  # Adjust this for placement
                #                     bbox_transform=ax.transAxes)  # 6x zoom
                # axins.plot(x1, y1, color='C3', alpha=0.6)  # Same data as main plot
                # axins.plot(x1, y1_scanfit, '--', color='r')
                # axins.fill_between(x1, lower_bound1, upper_bound1, facecolor='C3', alpha=0.4)
                # axins.plot(x2, y2, color='C0', alpha=0.6)  # Same data as main plot
                # axins.plot(x2, y2_scanfit, '--', color='b')
                # axins.fill_between(x2, lower_bound2, upper_bound2, facecolor='C0', alpha=0.4)

                # # Define zoomed-in region
                # x1, x2, y1, y2 = max(x2)-0.01, max(x2), max(y2)-0.03, max(y1)+0.03
                # axins.set_xlim(x1, x2)
                # axins.set_ylim(y1, y2)

                # # Customize inset tick labels
                # axins.yaxis.get_major_locator().set_params(nbins=7)
                # axins.ticklabel_format(useOffset=False, style='plain')
                # axins.xaxis.get_major_locator().set_params(nbins=3)
                # axins.tick_params(labelleft=True, labelbottom=True, labelsize=12)

                # # Mark the zoomed-in area on the main plot
                # mark_inset(ax, axins, loc1=1, loc2=2, fc="none", ec="black", linestyle="dashed")

        # ax.set_title(f'$\chi^2/\\text{{dof}}$={chi2_value1:.1f}/{dof1}', fontsize=25)
        
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(False)
        ax.legend(loc='best', fontsize=25)
        plt.savefig(os.path.join(plots, f'{date}', f'{dtype}_vs_time_{date}_run{run}.png'))
        # plt.show()

    def allan_deviation(self, timestamp, wavelength, run, name, xlabel, ylabel, date):
        """
        Calculate the Allan deviation of the measured wavelength/frequency.
        :param timestamp: Time data from the wavelengthmeter
        :param wavelength: Wavelength data from the wavelengthmeter
        :param run: Current run index
        :param name: Name of the measured quantity
        :param ave_unit: Average unit
        :param std_unit: Standard deviation unit
        :param xlabel: X-axis label
        :param ylabel: Y-axis label
        :param date: Date of the measurement
        """
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        x, wl, nu = [], [], []

        for i in self.number_of_runs(run):
            # Filter the data and convert wavelength to frequency
            x, wl = self.analyzer.filter_data(timestamp[i], wavelength[i]*1e-9)  # [s]
            nu = np.array([self.consts.c / wl[j] for j in range(len(wl))], dtype=np.float64)  # [Hz]
            nu_fit, slope, intercept, nu_weightedmean, residuals, nu_residualstd = self.analyzer.drift_fit(x, nu, 100)

            y = (nu - nu_weightedmean) / nu_weightedmean

            sample_interval = np.float64(x[-1] / len(y))  # Ensure high precision

            taus, adev, errors, _ = allantools.oadev(y, data_type="freq", rate=1/sample_interval)
            print("Optimal averaging time: ", taus[np.argmin(adev)])

            # Plot the Allan deviation
            ax.errorbar(taus, adev, yerr=errors, fmt="o", linestyle="-", color="C0", label="Allan Deviation", capsize=5)

        # Set logarithmic scale for both axes
        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel(xlabel, fontsize=25)
        ax.set_ylabel(ylabel, fontsize=25)
        ax.tick_params(axis='x', labelsize=25)
        ax.tick_params(axis='y', labelsize=25)
        # ax.set_ylim(1e5, 1e7)
        # ax.set_yticks([])
        ax.legend(loc='best', fontsize=25)
        ax.grid(which="both", linestyle="--")
        plt.savefig(os.path.join(plots, f'{date}', f'{name}_allan_deviation_mdev_{date}_run{run}.png'))
        # plt.show()
        
    def wavelength_frequency(self, wavelengthmeter_path, date, run, dtype):
        # Read wavelength data from Bristol wavelength meter
        timestamp, wavelength = self.reader.read_bristol(wavelengthmeter_path)

        if dtype == 'wavelength':
            name = r'\lambda'
            ave_unit = r'nm'
            std_unit = r'pm'
            self.plot_process(timestamp, wavelength, dtype, run, name, ave_unit, std_unit, 
                            r'Time (min)', r'Wavelength (nm)', date)
        elif dtype == 'frequency':
            name = r'\nu'
            ave_unit = r'GHz'
            std_unit = r'MHz'
            self.plot_process(timestamp, wavelength, dtype, run, name, ave_unit, std_unit,  
                            r'Time (min)', r'Frequency (GHz)', date)
            
        elif dtype == 'fft':
            xlabel = r'Frequency (Hz)'
            ylabel = r'Amplitude'
            file_name = f'FFT_{date}_run{run}.png'
            self.plot_fft(timestamp, wavelength, run, xlabel, ylabel, file_name)

        elif dtype == 'allan':
            name = r'Frequency'
            xlabel = r'Averaging Time, $\tau$ (s)'
            ylabel = r'Allan Deviation, $\sigma_\nu(\tau)$'
            self.allan_deviation(timestamp, wavelength, run, name, xlabel, ylabel, date)

if __name__ == "__main__":
    dir_path = os.path.join(
    os.path.expanduser('~'),  # Directory path on personal computer
    # 'D:',
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
    date_input = '02-26-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    wavelengthmeter_path = glob.glob(os.path.join(wavelengthmeter, date, '*.csv'))
    for i in range(1, 10):
        plotter.wavelength_frequency(wavelengthmeter_path, date, i, 'frequency')
        # plotter.wavelength_frequency(wavelengthmeter_path, date, i, 'allan')