import os, glob
import datetime as dt
from datetime import datetime
import numpy as np
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
        self.l = (7.5 - 0.159 * 2) * 1e-2  # Optical path length in [m]

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run)
    
    def process_physics(self, lambda_path, lockin_path, dtype, n, run):
        """
        Calculate ellipticity and Faraday rotation from measured data.
        :param lambda_path: Path to wavelength data
        :param lockin_path: Path to lock-in data
        :param dtype: Data type ('X' for in-phase, 'R' for magnitude)
        :param n: Skipping data points equal to n x Time Constant
        :param run: Current run index
        :return: Processed wavelength, detuning, ellipticity, and rotation angle
        """
        # Read wavelength data from Bristol wavelength meter
        B_t, Lambda = self.reader.read_bristol(lambda_path)

        # Process data based on the selected dtype
        if dtype == 'X':
            # Extract in-phase components and calculate ellipticity/angle
            para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f = self.reader.read_lockins(lockin_path)
            epsilon, epsilon_approx = self.analyzer.ellipticity(lockin_path, X1f, Xdc)
            theta = self.analyzer.angle(lockin_path, X1f, X2f, Xdc)
            m2f_theta = self.analyzer.modulated_angle(lockin_path, X1f, Xdc, Xm2f)
        elif dtype == 'R':
            # Extract magnitude components and calculate ellipticity/angle
            para, lockins_t, R1f, R2f, Rdc, Rm2f = self.analyzer.R_lockins(lockin_path)
            epsilon, epsilon_approx = self.analyzer.ellipticity(lockin_path, R1f, Rdc)
            theta = self.analyzer.angle(lockin_path, R1f, R2f, Rdc)
            m2f_theta = self.analyzer.modulated_angle(lockin_path, R1f, Rdc, Rm2f)

        # Initialize lists for processed data
        timestamp, wavelength, detuning, ellipticity, angle, m2f_angle = [], [], [], [], [], []

        # Process data for each run
        for i in self.number_of_runs(run):
            # Linear fit for laser drift
            frequency = [self.consts.c / Lambda[j] for j in range(len(Lambda))]  # [GHz]
            y_fit, slope, intercept, y_weightedmean, residuals, y_residualstd = self.analyzer.drift_fit(B_t[i], frequency[i], 100)

            # Filter and trim data for the current run
            B_t[i], Lambda[i] = self.analyzer.filter_data(B_t[i], Lambda[i])
            B_t[i], Lambda[i], lockins_t[i], epsilon_trimmed = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], epsilon[i])
            B_t[i], Lambda[i], lockins_t[i], theta_trimmed = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], theta[i])
            B_t[i], Lambda[i], lockins_t[i], m2f_theta_trimmed = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], m2f_theta[i])

            # Calculate intervals and averages
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[i], lockins_t[i], para[i][2], n)
            Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon_trimmed[l_idx])
            Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta_trimmed[l_idx])
            Lambd, m2f_th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], m2f_theta_trimmed[l_idx])

            # Append processed data to the respective lists
            timestamp.append(B_t[i][b_idx])  # [s]
            wavelength.append(Lambd)  # [m]
            ellipticity.append(ep*1e3)  # [mrad]
            angle.append(th*1e3)  # [mrad]
            m2f_angle.append(m2f_th*1e6)  # [μrad]
            detuning = y_weightedmean * 1e-6 - self.consts.K39_D2_Hz * 1e-6  # [MHz]
            detuning_std = y_residualstd * 1e-6  # [MHz]
            detuning_sem = detuning_std / np.sqrt(len(Lambd)) # [MHz]
            
        index = 1 if len(self.number_of_runs(run)) > 1 else 0
        # Initialize lists to store absorbance difference and refractive index difference
        alpha_diff_vapor, n_diff_vapor = [], []
        # Calculate absorbance and refractive index differences for vapor material
        alpha_diff_vapor.append(
            self.analyzer.absorbance_difference(ellipticity[index][1:], self.l)
            )  # Convert ellipticity to mrad for absorbance calculation
        n_diff_vapor.append(
            self.analyzer.refractive_indices_difference(angle[index][1:] * 1e6, self.l, wavelength[index])
            )  # Convert Faraday rotation to nanoradians for refractive index calculation

        return timestamp, wavelength, detuning, detuning_sem, ellipticity, angle, m2f_angle, alpha_diff_vapor, n_diff_vapor, index
    
    def raw_plot(self, lambda_path, lockin_path, dtype, n, run, temp, power, phytype, material, date):
        """
        Plot raw ellipticities and Faraday rotation angles measured by the lock-in amplifiers.
        :param lambda_path: Path to wavelength data
        :param lockin_path: Path to lock-in data
        :param dtype: Data type ('X' or 'R')
        :param n: Number of data points to average
        :param run: Current run index
        :param B: Magnetic field strength [G]
        :param power: Laser power [μW]
        :param phytype: Physical quantity type ('CD' or 'CB')
        :param material: Measurement material ('air', 'empty', 'vapor', etc.)
        """
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        # Import the processed data
        timestamp, wavelength, detuning, detuning_sem, ellipticity, angle, m2f_angle, alpha_diff_vapor, n_diff_vapor, index = \
            self.process_physics(lambda_path, lockin_path, dtype, n, run)
        
        # Calculate the average magnetic field and its variation
        B_mean, T_mean, B_std, T_std, _, _, _, _ = self.Bfield_and_temperature(gaussmeter_path, run)

        plot_params = {
            ('CD', 'vapor'): (timestamp[index] / 3600, ellipticity[index], r'$\epsilon_\text{vapor cell}$'),
            ('CB', 'vapor'): (timestamp[index] / 3600, angle[index], r'$\theta_\text{vapor cell}$'),
            ('modCB', 'vapor'): (timestamp[index] / 3600, m2f_angle[index], r'$\Delta\theta$'),
            # ('modCB', 'vapor'): (timestamp[index] / 3600, m2f_angle[index], r'$\Delta\theta/\Delta B_z$'),
            ('absorbance', 'vapor'): (timestamp[index][1:] / 60, alpha_diff_vapor[index], r'$\alpha_--\alpha_+$'),
            ('refractive', 'vapor'): (timestamp[index][1:] / 60, n_diff_vapor[index], r'$n_--n_+$'),
        }
        # Check if the combination of phytype and material exists in the mapping
        key = (phytype, material)
        if key in plot_params:
            # Retrieve plotting data (x-axis, y-axis, label) and plot on the axes
            x, y, label = plot_params[key]

            # calculate Δθ/ΔB
            # if key == ('modCB', 'vapor'):
                # y =  np.array(y) / (2 * B_std * 1e3) # [nrad/mG]

        # Fitting the data to a linear model
        y_fit, slope, intercept, y_weightedmean, residual, y_residualstd = self.analyzer.drift_fit(x, y, int(len(y)*0.1))
        y_linearfit, fitted_k, fitted_b, sigma_k, sigma_b, residuals, std, chi2_value, dof = self.analyzer.linear_fit(x, y, slope, intercept, y_residualstd)
        sem = std / np.sqrt(len(y))

        # the 1 sigma upper and lower analytic population bounds
        lower_bound = y_linearfit - std
        upper_bound = y_linearfit + std
        
        # plot Δθ vs time
        ax.plot(x, y, color='C0', alpha=1, lw=2, marker='o', markersize=5, label=f'$\\overline{{\\Delta\\theta}}$={round(y_weightedmean,2):.2f} μrad ± {round(sem*1e3)} nrad')
        ax.plot(x, y_linearfit, '--', color='b', lw=2, 
            label=f'$\\dot{{\\Delta\\theta}}$={round(fitted_k,2):.2f} μrad/h')
        ax.fill_between(x, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(std,2):.2f} μrad')

        # plot Δθ/ΔB vs time
        # ax.plot(x, y, color='C0', lw=2, marker='o', markersize=5, label=f'$\\overline{{\\Delta\\theta/\Delta B_z}}$={round(y_mean)} nrad/mG')
        # ax.plot(x, y_fit, '--', color='b', lw=2, 
        #     label=fr'$\nabla_t \frac{{\overline{{\Delta\theta}}}}{{\overline{{\Delta B_z}}}} = {round(slope,2):.2f} \,\mathrm{{nrad/(mG\cdot min)}}$')
        # ax.fill_between(x, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(y_std)} nrad/mG')
        
        # here we use the where argument to only fill the region where the
        # walker is above the population 1 sigma boundary
        ax.fill_between(x, upper_bound, y, where=y > upper_bound, fc='r', alpha=0.4, interpolate=True)
        ax.fill_between(x, lower_bound, y, where=y < lower_bound, fc='r', alpha=0.4, interpolate=True)

        B_ave, B_spread = self.B_field()
        
        # ax.set_title(f'$\chi^2/\\text{{dof}}$={chi2_value:.1f}/{dof}', fontsize=25)
        self.plot_settings(run, B_ave, B_spread, detuning, detuning_sem, temp, power, date, dtype, phytype)

    def two_axes_plot(self, lambda_path, lockin_path, dtype, n, run, B, power, phytype, material, date):
        """
        Plot background-subtracted ellipticities and optical rotation angles.
        :param lambda_path: Path to wavelength data
        :param lockin_path: Path to lock-in data
        :param dtype: Data type ('X' or 'R')
        :param n: Skipping data points equal to n x Time Constant
        :param run: Current run index
        :param B: Magnetic field strength [G]
        :param power: Laser power [μW]
        :param phytype: Physical quantity type ('CD', 'CB', 'absorbance', 'refractive index')
        :param material: Measurement material ('air', 'empty', 'vapor', etc.)
        """
        # Create a figure and axes for plotting
        fig, ax1 = plt.subplots(1, 1, figsize=(25.60, 14.40))
        # Create second Y-axis
        ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis

        timestamp, wavelength, detuning, detuning_std, ellipticity, angle, m2f_angle, alpha_diff_vapor, n_diff_vapor, index = \
            self.process_physics(lambda_path, lockin_path, dtype, n, run)
        
        plot_params = {
            ('CD', 'vapor'): (timestamp[index] / 3600, ellipticity[index], angle[index], \
                              r'$\Longleftarrow$$\epsilon_\text{vapor cell}$', \
                                r'$\theta_\text{vapor cell}$$\Longrightarrow$', \
                                  r'Ellipticity (mrad)', r'Faraday Rotation (mrad)', \
                                    f'[{dtype}]Ellipticity_and_Rotation_vapor_{date}_run{run}.png'),
            ('absorbance', 'vapor'): (timestamp[index][1:] / 3600, alpha_diff_vapor[index], n_diff_vapor[index], \
                                      r'$\Longleftarrow$$\alpha_--\alpha_+$', r'$n_--n_+$$\Longrightarrow$', \
                                        r'Absorbance difference (1/m)', r'Refractive indices difference ($\times10^{-6}$)', \
                                          f'[{dtype}]Absorbance_and_refractive_index_vapor_{date}_run{run}.png'),
        }
        B_ave, B_spread = self.B_field()

        # Check if the combination of phytype and material exists in the mapping
        key = (phytype, material)
        if key in plot_params:
            # Retrieve plotting data (x-axis, y-axis, label) and plot on the axes
            x, y1, y2, label_y1, label_y2, ylabel_y1, ylabel_y2, file_name = plot_params[key]

            y1_fit, y1_slope, y1_intercept, y1_weightedmean, residual, y1_residualstd = self.analyzer.drift_fit(x, y1, int(len(y1)*0.1))
            y2_fit, y2_slope, y2_intercept, y2_weightedmean, residual, y2_residualstd = self.analyzer.drift_fit(x, y2, int(len(y2)*0.1))

            y1_linearfit, fitted_k1, fitted_b1, sigma_k1, sigma_b1, y1_residuals, y1_std, chi2_value1, dof1 = self.analyzer.linear_fit(x, y1, y1_slope, y1_intercept, y1_residualstd)
            y2_linearfit, fitted_k2, fitted_b2, sigma_k2, sigma_b2, y2_residuals, y2_std, chi2_value2, dof2 = self.analyzer.linear_fit(x, y2, y2_slope, y2_intercept, y2_residualstd)

            y1_sem = y1_residualstd / np.sqrt(len(y1))
            y2_sem = y2_residualstd / np.sqrt(len(y2))

            y1_lower_bound = y1_linearfit - y1_std
            y1_upper_bound = y1_linearfit + y1_std

            y2_lower_bound = y2_linearfit - y2_std
            y2_upper_bound = y2_linearfit + y2_std

            ax1.plot(x, y1, color='C3', linestyle='-', linewidth=1, marker='o', markersize=5, \
                    label=f'$\\overline{{\\epsilon}}$={round(y1_weightedmean,2):.2f} mrad ± {round(y1_sem*1e3,2):.2f} μrad')
            ax1.plot(x, y1_fit, '--', color='r', lw=2, \
                    label=f'$\\dot\\epsilon$={round(fitted_k1*1e3,2):.2f} μrad/h')
            ax1.fill_between(x, y1_lower_bound, y1_upper_bound, facecolor='C3', alpha=0.4, label=f'$\\sigma_\epsilon$={round(y1_std*1e3,2):.2f} μrad')
            ax1.set_xlabel(r'Time (h)', fontsize=25)
            ax1.set_ylabel(ylabel_y1, fontsize=25, color='C3')
            ax1.tick_params(axis='x', labelsize=25)
            ax1.tick_params(axis='y', labelsize=25)
            # ax1.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.2f'))

            ax2.plot(x, y2, color='C0', linestyle='-', linewidth=1, marker='o', markersize=5, \
                    label=f'$\\overline{{\\theta}}$={round(y2_weightedmean,2):.2f} mrad ± {round(y2_sem*1e3,2):.2f} μrad')
            ax2.plot(x, y2_fit, '--', color='b', lw=2, 
                    label=f'$\\dot\\theta$={round(fitted_k2*1e3,2):.2f} μrad/h')
            ax2.fill_between(x, y2_lower_bound, y2_upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma_\\theta$={round(y2_std*1e3,2):.2f} μrad')
            ax2.set_ylabel(ylabel_y2, fontsize=25, color='C0')
            ax2.tick_params(axis='y', labelsize=25)
            # ax2.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.2f'))

            if datetime.strptime(date, "%m-%d-%Y") < reference_date:
                plt.title(f'$B_z$={B} G, $P$={power} μW @{date}', fontsize=25)
            else:
                B_mean, T_mean, B_std, T_std, _, _, _, _ = self.Bfield_and_temperature(gaussmeter_path, run)
                # ax1.set_title(f'$\chi_\\epsilon^2/\\text{{dof}}$={chi2_value1:.1f}/{dof1}', fontsize=25, pad=15, loc='left')
                # ax2.set_title(f'$\chi_\\theta^2/\\text{{dof}}$={chi2_value2:.1f}/{dof2}', fontsize=25, pad=15, loc='right')

            # plt.title(fr'$B_z$={B_ave:.3f}$\pm${B_spread:.3f} G, $\Delta$={frequency_shift} MHz, $T$={temp:.2f}°C, $P$={power} μW', fontsize=25)
        
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="best", fontsize=25)
        # plt.grid()
        save_path = os.path.join(plots, date, file_name)
        plt.savefig(save_path)
        # plt.show()

    def plot_settings(self, run, B_ave, B_variation, detuning, detuning_sem, temp, power, date, dtype, phytype):
        """
        Configure plot settings for consistent visualization.
        :param run: Current run index
        :param B: Magnetic field strength [G]
        :param power: Laser power [μW]
        :param dtype: Data type ('X' or 'R')
        :param phytype: Physical quantity type ('CD', 'CB', etc.)
        """
        plt.xlabel(r'Time (h)', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.grid(False)
        plt.legend(loc='best', fontsize=25)
        plot_labels = {
        'CD': r'Ellipticity (mrad)',
        'CB': r'Faraday Rotation (mrad)',
        'modCB': r'$\Delta\theta$ (μrad)',
        # 'modCB': r'$\Delta\theta/\Delta B_z$ (nrad/mG)',
        'absorbance': r'$\alpha_--\alpha_+$ (rad/m)',
        'refractive': r'$n_--n_+$ ($\times10^{-6}$)',
        }
        B_mean, T_mean, B_std, T_std, _, _, _, _ = self.Bfield_and_temperature(gaussmeter_path, run)
        plot_titles = f'$B_z$={round(-B_mean,3):.3f}±{round(B_std,3):.3f} G, $T$={round(T_mean,2):.2f}°C, $\\Delta$={round(detuning,2):.2f}±{round(detuning_sem,2)} MHz, $P$={power} μW'
        
        # plot_titles = fr'$B_z$={B_ave:.3f}$\pm${B_variation:.3f} G, $\Delta$={nu_shift} MHz, $T$={temp:.2f}°C, $P$={power} μW'

        file_names = {
            'CD': f'[{dtype}]Ellipticity_vs_Time_{date}_run{run}.png',
            'CB': f'[{dtype}]FR_vs_Time_{date}_run{run}.png',
            'modCB': f'[{dtype}]Mod_FR_vs_Time_{date}_run{run}.png',
            # 'modCB': f'[{dtype}]Mod_FR_per_mG_vs_Time_{date}_run{run}.png',
            'absorbance': f'[{dtype}]Absorbance_vs_Time_{date}_run{run}.png',
            'refractive': f'[{dtype}]Refractive_index_vs_Time_{date}_run{run}.png',
        }

        if phytype in plot_labels:
            plt.ylabel(plot_labels[phytype], fontsize=25)
            # plt.title(plot_titles, fontsize=25)
            file_name = file_names[phytype].format(dtype=dtype, date=date, run=run)
            save_path = os.path.join(plots, date, file_name)
            plt.savefig(save_path)

        print("Average B field: ", round(-B_mean,3))
        print("Average temperature: ", round(T_mean,2))
        print("Detuning: ", round(detuning,2))
        print("Detuning uncertainty: ", round(detuning_sem,2))
        print("Laser power: ", power)

        # plt.show()

    def Bfield_and_temperature(self, gaussmeter_path, run):
        timestamps, B0s, temps = self.reader.read_gaussmeter(gaussmeter_path)
        B_linearfit, B_slope, B_intercept, B_weightedmean, B_residuals, B_residualstd = self.analyzer.drift_fit(timestamps[run-1], B0s[run-1], 20)
        T_mean = np.mean(temps[run-1])
        T_std = np.std(temps[run-1], ddof=1) / np.sqrt(len(temps[run-1]))

        return B_weightedmean, T_mean, B_residualstd, T_std, B_linearfit, B_slope, B_intercept, B_residuals

    def B_field(self):
        B_max = np.array([5.2934, 5.2915, 5.2932, 5.2927, 5.2935])
        B_min = np.array([5.2848, 5.2834, 5.2842, 5.2839, 5.2837])

        # Compute the average field
        B_avg = np.round(0.5 * (np.mean(B_max) + np.mean(B_min)),3)

        # Compute mean absolute deviation (MAD)
        B_spread = np.round(0.5 * (np.mean(np.abs(B_max - B_avg)) + np.mean(np.abs(B_min - B_avg))),3)

        return B_avg, B_spread

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
    gaussmeter = os.path.join(K_vapor, 'Gaussmeter_data')
    lockins = os.path.join(K_vapor, 'Lockins_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data')
    # Define the reference date
    reference_date = datetime.strptime("02-09-2025", "%m-%d-%Y")

    plotter = Plot()
    date_input = '02-26-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    wavelengthmeter_path = glob.glob(os.path.join(wavelengthmeter, date, '*.csv'))
    gaussmeter_path = glob.glob(os.path.join(gaussmeter, date, '*.csv'))
    lockins_path = glob.glob(os.path.join(lockins, date, '*.lvm'))
    # for i in range(1,6):
    plotter.raw_plot(wavelengthmeter_path, lockins_path, 'X', 8, 4, 22.75, 200, 'modCB', 'vapor', date)
    plotter.two_axes_plot(wavelengthmeter_path, lockins_path, 'X', 6, 4, 22.75, 200, 'CD', 'vapor', date)

    FR_file = f'FaradayRotation_{date_input}.csv'
    # plotter.write(Bristol_path, Lockins_path, processed_path, FR_file, 'X', 5, 3, 22.00, 0.005, 41.0)