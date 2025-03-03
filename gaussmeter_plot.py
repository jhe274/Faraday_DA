import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze

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
        T_sem = np.std(temps[run], ddof=1) / np.sqrt(len(temps[run]))

        return T_mean, T_sem
    
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

    def plot_process(self, x, y1, y2, run, xlabel, ylabel_y1, ylabel_y2, file_name):
        # Create a figure and axes for plotting
        fig, ax1 = plt.subplots(1, 1, figsize=(25.60, 14.40))
        # Create second Y-axis
        ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis
        instrument_uncertainty = 0.0005
        instrument_resolution = 0.02*1e-3
        
        for i in self.number_of_runs(run):
            if i == run-1:
                time_factor = 60  # Convert time to min
                y1[i] = -y1[i]  # Invert the magnetic field data
                T_mean, T_sem = self.process_temperature(y2, i)

                # Fitting constant magnetic field measurements
                y1_linearfit, slope, intercept, y1_weightedmean, linear_residuals, y1_residualstd = self.analyzer.drift_fit(x[i], y1[i], 20)

                # Calculate the uncertainties
                sem = y1_residualstd / np.sqrt(len(y1[i]))
                sigma_mean = max(y1_weightedmean * instrument_uncertainty, sem)
                
                # Fitting constant magnetic field measurements with initial guesses
                # k0 = slope  
                # b0 = intercept
                # y_linearfit, fitted_k, fitted_b, sigma_k, sigma_b, residuals, std, chi2_value, dof = self.analyzer.linear_fit(x[i], y1[i], k0, b0, sigma_mean)
                # sigma_drift = max(sigma_k, y1_weightedmean * instrument_uncertainty / max(x[i]))
            
                # Apply band-pass filter around detected frequency
                # lowcut, highcut = dominant_freq - 0.1, dominant_freq + 0.1
                # y1_filtered = self.bandpass_filter(y1[i], lowcut, highcut, fs=10)

                # Fitting modulated magnetic field measurements
                freqs, fft_values, amplitude_spectrum, dominant_freq = self.analyzer.fft_peak(x[i], y1[i])
                k0 = slope # slope guess
                b0 = intercept # intercept guess
                a0 = 45 * 1e-3 # amplitude guess
                phi0 = 0 # phase guess  
                y_sinefit, fitted_k, fitted_b, fitted_a, fitted_phi, sigma_k, sigma_b, sigma_a, sigma_phi, residuals, std, chi2_value, dof = self.analyzer.drift_sine_fit(x[i], y1[i], k0, b0, a0, phi0, sigma_mean)
                sigma_drift = max(sigma_k, y1_weightedmean * instrument_uncertainty / max(x[i]))

                # Manual fitting of modulated magnetic field measurements
                # y_test = self.analyzer.drift_sine_wave(x[i], fitted_k, fitted_b, 80*1e-3, fitted_phi+np.pi)

                # Calculate the error bounds
                # lower_bound = y_linearfit - y1_residualstd
                # upper_bound = y_linearfit + y1_residualstd
                lower_bound = y_sinefit - y1_residualstd
                upper_bound = y_sinefit + y1_residualstd

                # Plot the magnetic field measurements
                ax1.plot(x[i]/time_factor, y1[i], color='C0', alpha=0.4, label=fr'$\dot{{B_z}}$={round(fitted_k*1e3*time_factor*60,1):.1f} mG/h')

                # plot the linear fit
                # ax1.plot(x[i]/time_factor, y_linearfit, '--', color='b', label=fr'$B_z$={round(y1_weightedmean,3):.3f} G ± {round(sigma_mean*1e3,1):.1f} mG')
                
                # Plot the modualted magnetic field fit
                ax1.plot(x[i]/time_factor, y_sinefit, '--', color='b', label=fr'$B_z$={round(y1_weightedmean,3):.3f} G ± {round(np.abs(fitted_a)*1e3,1):.1f} mG')

                # Plot the error bars
                ax1.fill_between(x[i]/time_factor, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(std*1e3,1):.1f} mG')
                # ax1.fill_between(x[i]/time_factor, upper_bound, y1[i], where=y1[i] > upper_bound, fc='red', alpha=0.4, interpolate=True)
                # ax1.fill_between(x[i]/time_factor, lower_bound, y1[i], where=y1[i] < lower_bound, fc='red', alpha=0.4, interpolate=True)

                ax1.set_xlabel(xlabel, fontsize=25)
                ax1.set_ylabel(ylabel_y1, fontsize=25)
                ax1.tick_params(axis='x', labelsize=25)
                ax1.tick_params(axis='y', labelsize=25)
                # ax1.set_title(f'$\chi^2/\\text{{dof}}$={chi2_value:.1f}/{dof}', fontsize=25)
                # Plot the temperature measurements
                ax2.plot(x[i]/time_factor, y2[i], label=f'$\\overline{{T}}$={round(T_mean,2):.2f}°C', color='black', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=2000)

                ax2.set_ylabel(ylabel_y2, fontsize=25)
                # ax2.set_yticks(np.arange(22.06, 22.07, 0.002))
                ax2.tick_params(axis='y', labelsize=25)

                # --------- Add Inset Zoomed Plot ---------
                axins = inset_axes(ax1, width="50%", height="50%", 
                                   bbox_to_anchor=(-0.18, -0.2, 0.5, 0.5),  # Adjust this for placement
                                    bbox_transform=ax1.transAxes)  # 6x zoom
                axins.plot(x[i]/time_factor, y1[i], color='C0', alpha=0.6)  # Same data as main plot
                axins.plot(x[i]/time_factor, y_sinefit, '--', color='b')
                axins.fill_between(x[i]/time_factor, lower_bound, upper_bound, facecolor='C0', alpha=0.4)

                # Define zoomed-in region
                x1, x2, y1, y2 = 0, 10/time_factor, fitted_b-y1_residualstd-0.1, fitted_b+y1_residualstd+0.1
                axins.set_xlim(x1, x2)
                axins.set_ylim(y1, y2)

                # Customize inset tick labels
                axins.yaxis.get_major_locator().set_params(nbins=7)
                axins.xaxis.get_major_locator().set_params(nbins=5)
                axins.tick_params(labelleft=True, labelbottom=True, labelsize=12)

                # Mark the zoomed-in area on the main plot
                mark_inset(ax1, axins, loc1=1, loc2=3, fc="none", ec="black", linestyle="dashed")
            else:
                time_factor = 60  # Convert time to hours
                y1[i] = -y1[i]  # Invert the magnetic field data
                T_mean, T_std = self.process_temperature(y2, i)

                # Fitting constant magnetic field measurements
                y1_linearfit, slope, intercept, y1_weightedmean, linear_residuals, y1_residualstd = self.analyzer.drift_fit(x[i], y1[i], 20)
                
                # Calculate the uncertainties
                sem = y1_residualstd / np.sqrt(len(y1[i]))
                sigma_mean = max(y1_weightedmean * instrument_uncertainty, sem)

                # Fitting constant magnetic field measurements with initial guesses
                k0 = slope
                b0 = intercept
                y_linearfit, fitted_k, fitted_b, sigma_k, sigma_b, chi2_value, dof = self.analyzer.linear_fit(x[i], y1[i], k0, b0, sigma_mean)
                sigma_drift = max(sigma_k, y1_weightedmean * instrument_uncertainty / max(x[i]))

                # Calculate the error bounds
                lower_bound = fitted_k * x[i] + fitted_b - y1_residualstd
                upper_bound = fitted_k * x[i] + fitted_b + y1_residualstd

                # Plot the magnetic field measurements
                ax1.plot(x[i]/time_factor, y1[i], color='C3', alpha=0.4, label=fr'$\dot{{B_z}}$={round(fitted_k*1e3*time_factor,2):.2f} mG/min')

                # plot the linear fit
                ax1.plot(x[i]/time_factor, y_linearfit, '--', color='r', label=fr'$B_z$={round(y1_weightedmean,3):.3f} G ± {round(sigma_mean*1e3,1):.1f} mG')

                # plot the error bars
                ax1.fill_between(x[i]/time_factor, lower_bound, upper_bound, facecolor='C3', alpha=0.4, label=f'$\\sigma$={round(y1_residualstd*1e3,1):.1f} mG')
                # ax1.fill_between(x[i]/time_factor, upper_bound, y1[i], where=y1[i] > upper_bound, fc='r', alpha=0.4, interpolate=True)
                # ax1.fill_between(x[i]/time_factor, lower_bound, y1[i], where=y1[i] < lower_bound, fc='r', alpha=0.4, interpolate=True)

                # Plot the temperature measurements
                ax2.plot(x[i]/time_factor, y2[i], label=f'$\\overline{{T}}$={round(T_mean,2):.2f}°C', color='C3', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=300)

        # plt.title(f'$\chi^2/\\text{{dof}}$={chi2_value:.1f}/{dof}', fontsize=25)
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=20)
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
    gaussmeter = os.path.join(K_vapor, 'Gaussmeter_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '03-02-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    gaussmeter_path = glob.glob(os.path.join(gaussmeter, date, '*.csv'))
    # for i in range(5,8):
    plotter.gaussmter_vs_time(gaussmeter_path, 15)