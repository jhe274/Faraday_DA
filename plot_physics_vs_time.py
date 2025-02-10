import os, glob
import datetime as dt
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
            y_fit, slope, intercept, y_mean, residuals, y_std = self.ld.drift_fit(B_t[i], frequency[i])

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
            detuning.append(self.consts.c / wavelength[i-run+1] * 1e-9 - self.consts.K39_D2_Hz * 1e-9)  # [GHz]
            ellipticity.append(ep*1e3)  # [mrad]
            angle.append(th*1e3)  # [mrad]
            m2f_angle.append(m2f_th*1e6)  # [μrad]

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

        return timestamp, wavelength, detuning, ellipticity, angle, m2f_angle, alpha_diff_vapor, n_diff_vapor, round(y_mean*1e-6-self.consts.K39_D2_Hz * 1e-6), index
    
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
        timestamp, wavelength, detuning, ellipticity, angle, m2f_angle, alpha_diff_vapor, n_diff_vapor, frequency_shift, index = \
            self.process_physics(lambda_path, lockin_path, dtype, n, run)
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        plot_params = {
            ('CD', 'vapor'): (timestamp[index] / 60, ellipticity[index], r'$\epsilon_\text{vapor cell}$'),
            ('CB', 'vapor'): (timestamp[index] / 60, angle[index], r'$\theta_\text{vapor cell}$'),
            ('modCB', 'vapor'): (timestamp[index] / 60, m2f_angle[index], r'$\theta_\text{m2f}$'),
            ('absorbance', 'vapor'): (timestamp[index][1:] / 60, alpha_diff_vapor[index], r'$\alpha_--\alpha_+$'),
            ('refractive', 'vapor'): (timestamp[index][1:] / 60, n_diff_vapor[index], r'$n_--n_+$'),
        }
        # Check if the combination of phytype and material exists in the mapping
        key = (phytype, material)
        if key in plot_params:
            # Retrieve plotting data (x-axis, y-axis, label) and plot on the axes
            x, y, label = plot_params[key]
        
        y_fit, slope, intercept, y_mean, residuals, y_std = self.analyzer.drift_fit(x, y)

        # the 1 sigma upper and lower analytic population bounds
        lower_bound = slope*x + intercept - y_std
        upper_bound = slope*x + intercept + y_std

        # ax.scatter(x, y, color='C0', s=70)
        ax.plot(x, y, lw=2, label=label)
        ax.plot(x, y_fit, '--', color='C0', lw=2, 
            label=r'Sample mean')
        ax.fill_between(x, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'1$\sigma$ range')
        
        # here we use the where argument to only fill the region where the
        # walker is above the population 1 sigma boundary
        ax.fill_between(x, upper_bound, y, where=y > upper_bound, fc='red', alpha=0.4, interpolate=True)
        ax.fill_between(x, lower_bound, y, where=y < lower_bound, fc='red', alpha=0.4, interpolate=True)

        B_ave, B_spread = self.B_field()
        print('Average measured change of polarization rotation:', y_mean)
        print('Standard deviation of residuals:', y_std)
        self.plot_settings(run, B_ave, B_spread, frequency_shift, temp, power, date, dtype, phytype)

    def two_axes_plot(self, lambda_path, lockin_path, dtype, n, run, temp, power, phytype, material, date):
        timestamp, wavelength, detuning, ellipticity, angle, m2f_angle, alpha_diff_vapor, n_diff_vapor, frequency_shift, index = \
            self.process_physics(lambda_path, lockin_path, dtype, n, run)
        # Create a figure and axes for plotting
        fig, ax1 = plt.subplots(1, 1, figsize=(25.60, 14.40))
        # Create second Y-axis
        ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis
        
        plot_params = {
            ('CD', 'vapor'): (timestamp[index] / 60, ellipticity[index], angle[index], \
                              r'$\Longleftarrow$$\epsilon_\text{vapor cell}$', \
                                r'$\theta_\text{vapor cell}$$\Longrightarrow$', \
                                  r'Ellipticity (mrad)', r'Faraday Rotation (mrad)', \
                                    f'[{dtype}]Ellipticity_and_Rotation_vapor_{date}_run{run}.png'),
            ('absorbance', 'vapor'): (timestamp[index][1:] / 60, alpha_diff_vapor[index], n_diff_vapor[index], \
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

            ax1.scatter(x, y1, color='r', label=label_y1, s=70)
            ax1.plot(x, y1, color='r', linestyle='-', linewidth=2)
            ax1.set_xlabel(r'Time (min)', fontsize=25)
            # ax1.set_xticks(np.arange(-5, 6, 1))
            ax1.set_ylabel(ylabel_y1, fontsize=25)
            ax1.tick_params(axis='x', labelsize=25)
            ax1.tick_params(axis='y', labelsize=25)
            ax1.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.0f'))

            ax2.scatter(x, y2, color='b', label=label_y2, s=70)
            ax2.plot(x, y2, '.', color='b', linestyle='-', linewidth=2)
            ax2.set_ylabel(ylabel_y2, fontsize=25)
            ax2.tick_params(axis='y', labelsize=25)
            ax2.get_yaxis().set_major_formatter(plt.FormatStrFormatter('%.1f'))

            plt.title(fr'$B_z$={B_ave:.3f}$\pm${B_spread:.3f} G, $\Delta$={frequency_shift} MHz, $T$={temp:.2f}°C, $P$={power} μW', fontsize=25)
        
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="lower left", fontsize=25)
        plt.grid(False)
        save_path = os.path.join(plots, date, file_name)
        plt.savefig(save_path)
        plt.show()

    def plot_settings(self, run, B_ave, B_variation, nu_shift, temp, power, date, dtype, phytype):
        """
        Configure plot settings for consistent visualization.
        :param run: Current run index
        :param B: Magnetic field strength [G]
        :param power: Laser power [μW]
        :param dtype: Data type ('X' or 'R')
        :param phytype: Physical quantity type ('CD', 'CB', etc.)
        """
        plt.xlabel(r'Time (min)', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.grid(False)
        plt.legend(loc='best', fontsize=25)
        plot_labels = {
        'CD': r'Ellipticity (mrad)',
        'CB': r'Faraday Rotation (mrad)',
        'modCB': r'$\Delta\theta$ (μrad)',
        'absorbance': r'$\alpha_--\alpha_+$ (rad/m)',
        'refractive': r'$n_--n_+$ ($\times10^{-6}$)',
        }
        plot_titles = fr'$B_z$={B_ave:.3f}$\pm${B_variation:.3f} G, $\Delta$={nu_shift} MHz, $T$={temp:.2f}°C, $P$={power} μW'
        file_names = {
            'CD': f'[{dtype}]Ellipticity_vs_Time_{date}_run{run}.png',
            'CB': f'[{dtype}]FR_vs_Time_{date}_run{run}.png',
            'modCB': f'[{dtype}]Mod_FR_vs_Time_{date}_run{run}.png',
            'absorbance': f'[{dtype}]Absorbance_vs_Time_{date}_run{run}.png',
            'refractive': f'[{dtype}]Refractive_index_vs_Time_{date}_run{run}.png',
        }

        if phytype in plot_labels:
            plt.ylabel(plot_labels[phytype], fontsize=25)
            plt.title(plot_titles, fontsize=25)
            
            file_name = file_names[phytype].format(dtype=dtype, date=date, run=run)
            save_path = os.path.join(plots, date, file_name)
            plt.savefig(save_path)

        plt.show()

    def write(self, lambda_path, lockin_path, folder_path, filename, dtype, n, run, T, B, P):
        """
        Write the processed data to a CSV file.

        :param lambda_path: Path to wavelength data
        :param lockin_path: Path to lock-in data
        :param folder_path: Path to the folder where the CSV file will be saved
        :param filename: Name of the output CSV file
        :param dtype: Data type ('X' or 'R')
        :param n: Skipping data points equal to n x Time Constant
        :param run: Current run index
        :param T: Temperature [°C]
        :param B: Longitudinal magnetic field [G]
        :param P: Laser power [μW]
        """
        # Extract processed data required for writing to CSV
        x0, x, CD_empty, CB_empty, CD_vapor, CB_vapor, CD_K, CB_K = \
            self.background_subtraction(lambda_path, lockin_path, dtype, n, run)
        
        # Data to be written (wavelength, ellipticity, Faraday rotation)
        data = [x0[0], CD_vapor, CB_vapor]

        try:
            # Counter for handling duplicate filenames
            counter = 1
            original_filename = filename

            # Ensure the folder for saving the file exists
            folder_path = os.path.join(folder_path, date_input)
            os.makedirs(folder_path, exist_ok=True)

            # Check if the file already exists and create a new unique filename if necessary
            while os.path.isfile(os.path.join(folder_path, filename)):
                filename = f"{original_filename.split('.')[0]}_{counter}.csv"
                counter += 1

            # Construct the full path for the output file
            file_path = os.path.join(folder_path, filename)

            # Write data to the CSV file
            with open(file_path, "w") as file:
                # Write metadata (date, temperature, field, and power) as the first lines
                for attribute, value in zip(
                    ['Date (MM-DD-YYYY)', 'Temperature (°C)', 'Longitudinal magnetic field (G)', 'Power (microW)'],
                    [date_input, T, B, P]
                ):
                    file.write(f'{attribute}, {value}\n')

                # Write the header for the data columns
                header = 'Wavelength (m), Ellipticity (radian), Faraday rotation (radian)\n'
                file.write(header)

                # Write the actual data rows
                for row in list(zip(*data)):
                    file.write(','.join(map(str, row)) + '\n')

        except Exception as e:
            # Handle any exceptions that occur during the file writing process
            print(f"An error occurred while saving data to the file: {e}")

    def B_field(self):
        B_max = np.array([5.2934, 5.2915, 5.2932, 5.2927, 5.2935])
        B_min = np.array([5.2848, 5.2834, 5.2842, 5.2839, 5.2837])

        # Compute the average field
        B_avg = np.round(0.5 * (np.mean(B_max) + np.mean(B_min)),3)
        print("Average Magnetic Field:", B_avg)

        # Compute mean absolute deviation (MAD)
        B_spread = np.round(0.5 * (np.mean(np.abs(B_max - B_avg)) + np.mean(np.abs(B_min - B_avg))),3)
        print("Average Magnetic Field Spread (Variation):", B_spread)

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
    lockins = os.path.join(K_vapor, 'Lockins_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data')
    
    plotter = Plot()
    date_input = '02-09-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    wavelengthmeter_path = glob.glob(os.path.join(wavelengthmeter, date, '*.csv'))
    lockins_path = glob.glob(os.path.join(lockins, date, '*.lvm'))
    # plotter.raw_plot(Bristol_path, Lockins_path, 'R', 8, 1, 22.75, 270, 'modCB', 'vapor', date)
    plotter.two_axes_plot(wavelengthmeter_path, lockins_path, 'X', 5, 1, 22.75, 395, 'CD', 'vapor', date)


    FR_file = f'FaradayRotation_{date_input}.csv'
    # plotter.write(Bristol_path, Lockins_path, processed_path, FR_file, 'X', 5, 3, 22.00, 0.005, 41.0)