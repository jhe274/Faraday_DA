import os, glob
import datetime as dt
import numpy as np
import scipy.special
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
        self.l = (7.5 - 0.159 * 2) * 1e-2  # Optical path length in meters

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run+1)
    
    def ellipticity_and_angle(self, lambda_path, lockin_path, dtype, n, run):
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
            para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xmod, Ymod = self.reader.read_lockins(lockin_path)
            epsilon, epsilon_approx = self.analyzer.ellipticity(lockin_path, X1f, Xdc)
            theta = self.analyzer.angle(lockin_path, X1f, X2f, Xdc)
        elif dtype == 'R':
            # Extract magnitude components and calculate ellipticity/angle
            para, lockins_t, R1f, R2f, Rdc = self.analyzer.R_lockins(lockin_path)
            epsilon, epsilon_approx = self.analyzer.ellipticity(lockin_path, R1f, Rdc)
            theta = self.analyzer.angle(lockin_path, R1f, R2f, Rdc)

        # Initialize lists for processed data
        wavelength, detuning, ellipticity, angle = [], [], [], []

        # Process data for each run
        for i in self.number_of_runs(run):
            # Filter and trim data for the current run
            B_t[i], Lambda[i] = self.analyzer.filter_data(B_t[i], Lambda[i])
            B_t[i], Lambda[i], lockins_t[i], epsilon_trimmed = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], epsilon[i])
            B_t[i], Lambda[i], lockins_t[i], theta_trimmed = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], theta[i])

            # Calculate intervals and averages for ellipticity and angle
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[i], lockins_t[i], para[i][2], n)
            Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon_trimmed[l_idx])
            Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta_trimmed[l_idx])

            # Append processed data to the respective lists
            wavelength.append(Lambd)  # [m]
            detuning.append(self.consts.c / wavelength[i-run+1] * 1e-9 - self.consts.K39_D2_Hz * 1e-9)  # [GHz]
            ellipticity.append(ep)  # [rad]
            angle.append(th)  # [rad]

        return wavelength, detuning, ellipticity, angle
    
    def background_subtraction(self, lambda_path, lockin_path, dtype, n, run):
        """
        Subtract background from ellipticity and Faraday rotation measurements.
        :param lambda_path: Path to wavelength data
        :param lockin_path: Path to lock-in data
        :param dtype: Data type ('X' or 'R')
        :param n: Skipping data points equal to n x Time Constant
        :param run: Current run index
        :return: Background-subtracted data
        """
        # Calculate ellipticity and angle for each run
        wavelength, detuning, ellipticity, angle = self.ellipticity_and_angle(lambda_path, lockin_path, dtype, n, run)

        # Initialize lists to store background-subtracted values
        CD_empty, CD_vapor, CD_K, CB_empty, CB_vapor, CB_K = [], [], [], [], [], []

        if len(self.number_of_runs(run)) > 2:
            # Subtract air/empty cell as background
            for idx, x_val in enumerate(wavelength[0]):
                idx_empty = np.argmin(np.abs(wavelength[1] - x_val))
                idx_vapor = np.argmin(np.abs(wavelength[2] - x_val))
                CD_empty.append(ellipticity[1][idx_empty] - ellipticity[0][idx])
                CD_vapor.append(ellipticity[2][idx_vapor] - angle[0][idx])
                CB_empty.append(angle[1][idx_empty] - angle[0][idx])
                CB_vapor.append(angle[2][idx_vapor] - angle[0][idx])
            for idx, x_val in enumerate(wavelength[1]):
                idx_K = np.argmin(np.abs(wavelength[2] - x_val))
                CD_K.append(ellipticity[2][idx_K] - ellipticity[1][idx])
                CB_K.append(angle[2][idx_K] - angle[1][idx])
        else:
            # Subtract air as background
            for idx, x_val in enumerate(wavelength[0]):
                idx_vapor = np.argmin(np.abs(wavelength[1] - x_val))
                CD_vapor.append(ellipticity[1][idx_vapor] - ellipticity[0][idx])
                CB_vapor.append(angle[1][idx_vapor] - angle[0][idx])

        return wavelength, detuning, CD_empty, CB_empty, CD_vapor, CB_vapor, CD_K, CB_K
    
    def raw_plot(self, lambda_path, lockin_path, dtype, n, run, B, power, phytype, material, date):
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
        wavelength, detuning, ellipticity, angle = self.ellipticity_and_angle(lambda_path, lockin_path, dtype, n, run)
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        # Plot data for each run based on the physical quantity type and material
        for i in self.number_of_runs(run):
            if phytype == 'CD':
                if material == 'air':
                    ax.plot(detuning[0], ellipticity[0][1:], label=r'$\epsilon_\text{air}$')
                elif material == 'empty' and len(self.number_of_runs(run)) > 2:
                    ax.plot(detuning[1], ellipticity[1][1:], label=r'$\epsilon_\text{empty cell}$')
                elif material == 'vapor':
                    index = 2 if len(self.number_of_runs(run)) > 2 else 1
                    ax.plot(detuning[index], ellipticity[index][1:], label=r'$\epsilon_\text{vapor cell}$')
            elif phytype == 'CB':
                if material == 'air':
                    ax.plot(detuning[0], angle[0][1:], label=r'$\theta_\text{air}$')
                elif material == 'empty' and len(self.number_of_runs(run)) > 2:
                    ax.plot(detuning[1], angle[1][1:], label=r'$\theta_\text{empty cell}$')
                elif material == 'vapor':
                    index = 2 if len(self.number_of_runs(run)) > 2 else 1
                    ax.plot(detuning[index], angle[index][1:], label=r'$\theta_\text{vapor cell}$')
        
        self.plot_settings(n, B, power, date, dtype, phytype)

    def background_subtracted_plot(self, lambda_path, lockin_path, dtype, n, run, B, power, phytype, material, date):
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
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        # Perform background subtraction and retrieve processed data
        wavelength, detuning, CD_empty, CB_empty, CD_vapor, CB_vapor, CD_K, CB_K = \
            self.background_subtraction(lambda_path, lockin_path, dtype, n, run)

        # Initialize lists to store absorbance difference and refractive index difference
        alpha_diff_vapor, n_diff_vapor = [], []

        # Calculate absorbance and refractive index differences for vapor material
        for i in range(len(detuning[0])):
            alpha_diff_vapor.append(
                self.analyzer.absorbance_difference(CD_vapor[i] * 1e3, self.l)
            )  # Convert ellipticity to mrad for absorbance calculation
            n_diff_vapor.append(
                self.analyzer.refractive_indices_difference(CB_vapor[i] * 1e9, self.l, wavelength[0][i])
            )  # Convert Faraday rotation to nanoradians for refractive index calculation

        # Define the mapping for plotting parameters based on physical type and material
        plot_params = {
            ('CD', 'empty'): (detuning[0], CD_empty, r'$\epsilon_\text{empty cell}-\epsilon_\text{air}$'),
            ('CD', 'vapor'): (detuning[0], CD_vapor, r'$\epsilon_\text{vapor cell}-\epsilon_\text{air}$'),
            ('CD', 'K'): (detuning[1], CD_K, r'$\epsilon_\text{vapor cell}-\epsilon_\text{empty cell}$'),
            ('CB', 'empty'): (detuning[0], CB_empty, r'$\theta_\text{empty cell}-\theta_\text{air}$'),
            ('CB', 'vapor'): (detuning[0], CB_vapor, r'$\theta_\text{vapor cell}-\theta_\text{air}$'),
            ('CB', 'K'): (detuning[1], CB_K, r'$\theta_\text{vapor cell}-\theta_\text{empty cell}$'),
            ('absorbance', 'vapor'): (detuning[0], alpha_diff_vapor, r'$\alpha_\text{vapor cell}-\alpha_\text{air}$'),
            ('refractive index', 'vapor'): (detuning[0], n_diff_vapor, r'$n_\text{vapor cell}-n_\text{air}$'),
        }

        # Check if the combination of phytype and material exists in the mapping
        key = (phytype, material)
        if key in plot_params:
            # Retrieve plotting data (x-axis, y-axis, label) and plot on the axes
            x, y, label = plot_params[key]
            ax.plot(x, y, '.', color='red', label=label, markersize=2)

        # Optionally plot peaks and valleys (commented out in the current implementation)
        # self.peaks_valleys_plot(x[0], CB_vapor)

        # Apply consistent plot settings (e.g., labels, title, grid) using a helper method
        self.plot_settings(run, B, power, date, dtype, phytype)

    def peaks_valleys_plot(self, x, y):
        """
        Find peaks and valleys in the ellipticity/FR
        """
        peaks, _ = find_peaks(np.array(y), prominence=20, height=(200, 600))
        valleys, _ = find_peaks(-1*np.array(y), prominence=10, height=(-600, -200))

        for idx, indices in enumerate([peaks, valleys]):
            for k in indices:
                plt.vlines(x=x[k], ymin=0, ymax=y[k], linestyle=':', color='purple')

    def curve_fitting(self):
        l = (7.5-0.159*2)*1e-2                                                                                                          # [m]
        nu_D1 = 389286.058716 * 1e9
        nu_D2 = 391016.17003 * 1e9

        def FR(nu, Kn, T, B, P, const):
            delta_nu_D2 = nu - nu_D2
            delta_nu_D1 = nu - nu_D1
            delta_doppler_D2 = self.theory.doppler_broad(nu_D2, T)
            delta_doppler_D1 = self.theory.doppler_broad(nu_D1, T)
            term1 = (
                (7*(delta_nu_D2**2 - delta_doppler_D2**2/4) / (delta_nu_D2**2 + delta_doppler_D2**2/4)**2) + 
                (4*(delta_nu_D1**2 - delta_doppler_D1**2/4) / (delta_nu_D1**2 + delta_doppler_D1**2/4)**2) - 
                (2*(delta_nu_D1*delta_nu_D2) / ((delta_nu_D2-delta_doppler_D2) * (delta_nu_D1-delta_doppler_D1))**2)) / (3 * self.consts.h)
            # term2 = np.sign(B) * ((nu / (nu_D1 * (nu - nu_D1))) - (nu / (nu_D2 * (nu - nu_D2)))) / (self.consts.k_B * T)
            
            diamagnetic_theta = self.consts.mu_B * B*1e-4 * (term1)
            paramagnetic_thetea = P * (
                    (delta_nu_D2 / ((delta_nu_D2 - delta_doppler_D2) ** 2)) -
                    (delta_nu_D1 / ((delta_nu_D1 - delta_doppler_D1) ** 2))
                )
            # paramagnetic_thetea = P * (
            #         (delta_nu_D2 / ((delta_nu_D2 - doppler_broad_D2)**2 + (doppler_broad_D2**2)/4)) -
            #         (delta_nu_D1 / ((delta_nu_D1 - doppler_broad_D1)**2 + (doppler_broad_D1**2)/4))
            #     )
            theta = self.consts.alpha * Kn*1e14 * l * (diamagnetic_theta + paramagnetic_thetea) * 1e6 + const                                                        # [microrad]
            return theta

        initial_guess = [1.47, 19.5, -5.103, -0.002, -60]
        bounds = ([.1, 15, -5.2, -1, -100], [5, 25, -5., 1, 100])
        # params, covariance = curve_fit(FR, x, y, p0=initial_guess, bounds=bounds)
        # print(params)
        # Kn, T, Bz, PK, const = np.round(params,3)
        # plt.plot(x*1e-9 - nu_D2 * 1e-9, FR(x, Kn, T, Bz, PK, const), '--', color='red', label='Curve fit')
        # plt.plot(x*1e-9 - nu_D2 * 1e-9, FR(x, 1.474, 19.497, -5.103, -.002, -60), '--', color='green', label='Manual fit')

    def plot_settings(self, run, B, power, date, dtype, phytype):
        """
        Configure plot settings for consistent visualization.
        :param run: Current run index
        :param B: Magnetic field strength [G]
        :param power: Laser power [μW]
        :param dtype: Data type ('X' or 'R')
        :param phytype: Physical quantity type ('CD', 'CB', etc.)
        """
        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        # plt.xticks(np.arange(-.8, 1.2, .1), fontsize=25)
        plt.yticks(fontsize=25)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(False)
        # plt.legend(loc='best', fontsize=25)
        plot_labels = {
        'CD': r'Ellipticity (μrad.)',
        'CB': r'Faraday Rotation (μrad.)',
        'absorbance': r'$\alpha_--\alpha_+$ (1/mm)',
        'refractive index': r'$n_--n_+$ ($\times10^{-9}$)',
        }

        plot_titles = {
            'CD': f'Ellipticity vs Frequency, $B_z$={B} G, $P$={power} μW @{date}',
            'CB': f'Faraday Rotation vs Frequency, $B_z$={B} G, $P$={power} μW @{date}',
            'absorbance': f'Absorbance vs Frequency, $B_z$={B} G, $P$={power} μW @{date}',
            'refractive index': f'Refractive index vs Frequency, $B_z$={B} G, $P$={power} μW @{date}',
        }

        file_names = {
            'CD': f'[{dtype}]Ellipticity_vs_Frequency_{date}_run{run}-{run+1}.pdf',
            'CB': f'[{dtype}]FR_vs_Frequency_{date}_run{run}-{run+1}.pdf',
            'absorbance': f'[{dtype}]Absorbance_vs_Frequency_{date}_run{run}-{run+1}.pdf',
            'refractive index': f'[{dtype}]Refractive_index_vs_Frequency_{date}_run{run}-{run+1}.pdf',
        }

        if phytype in plot_labels:
            plt.ylabel(plot_labels[phytype], fontsize=25)
            plt.title(plot_titles[phytype], fontsize=25)
            
            file_name = file_names[phytype].format(dtype=dtype, date=date, run=run)
            save_path = os.path.join(Plots, date, file_name)
            plt.savefig(save_path)

        # plt.title(rf'$n={Kn}\times10^{{14}}\text{{m}}^3$, $T={T}^\circ$C, $B_z={Bz}$G, $P=.2\%$, $\theta_\text{{offset}}={const}μ\text{{rad}}$', fontsize=25)
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
    K_vapor = os.path.join(dir_path, 'K_vapor_cell')
    Bristol = os.path.join(K_vapor, 'Bristol_data')
    Lockins = os.path.join(K_vapor, 'Lockins_data')
    Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data')
    
    plotter = Plot()
    date_input = '06-07-2024'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
    Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
    # plotter.raw_plot(Bristol_path, Lockins_path, 'X', 5, 11, -6.105, 0.5, 'CD', 'air')
    plotter.background_subtracted_plot(Bristol_path, Lockins_path, 'X', 5, 1, -6.05, 402.3, 'refractive index', 'vapor', date)

    FR_file = f'FaradayRotation_{date_input}.csv'
    # plotter.write(Bristol_path, Lockins_path, processed_path, FR_file, 'X', 5, 3, 22.00, 0.005, 41.0)