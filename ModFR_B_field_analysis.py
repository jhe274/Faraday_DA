import os, glob
import datetime as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from constants import Constants as Consts
from theory_calculations import Theory
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze

class FaradayRotation:
    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()
        self.l = (7.5 - 0.159 * 2) * 1e-2  # Optical path length in [m]

    def process_physics(self, lambda_path, lockin_path, dtype, n, run):
        """
        Calculate ellipticity and Faraday rotation from measured data.
        :param lambda_path: Path to wavelength data
        :param lockin_path: Path to lock-in data
        :param dtype: Data type ('X' for in-phase, 'R' for magnitude)
        :param n: Skipping data points equal to n x Time Constant
        :param run: Current run index
        :return: Processed wavelength, frequency, ellipticity, and rotation angle
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
        timestamp, wavelength, freq, ellipticity, angle, m2f_angle = [], [], [], [], [], []

        # Filter data for the current run
        B_t[run-1], Lambda[run-1] = self.analyzer.filter_data(B_t[run-1], Lambda[run-1])

        # Linear fit for laser drift
        freq = [self.consts.c / Lambda[run-1][j] for j in range(len(Lambda[run-1]))]  # [GHz]
        y_fit, slope, intercept, y_weightedmean, residuals, y_residualstd = self.analyzer.drift_fit(B_t[run-1], freq, 100)

        # Filter and trim data for the current run
        B_t[run-1], Lambda[run-1] = self.analyzer.filter_data(B_t[run-1], Lambda[run-1])
        B_t[run-1], Lambda[run-1], lockins_t[run-1], epsilon_trimmed = self.analyzer.trim_data(B_t[run-1], Lambda[run-1], lockins_t[run-1], epsilon[run-1])
        B_t[run-1], Lambda[run-1], lockins_t[run-1], theta_trimmed = self.analyzer.trim_data(B_t[run-1], Lambda[run-1], lockins_t[run-1], theta[run-1])
        B_t[run-1], Lambda[run-1], lockins_t[run-1], m2f_theta_trimmed = self.analyzer.trim_data(B_t[run-1], Lambda[run-1], lockins_t[run-1], m2f_theta[run-1])

        # Calculate intervals and averages
        l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[run-1], lockins_t[run-1], para[run-1][0], n)
        Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[run-1], Lambda[run-1][b_idx], epsilon_trimmed[l_idx])
        Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[run-1], Lambda[run-1][b_idx], theta_trimmed[l_idx])
        Lambd, m2f_th = self.analyzer.calculate_averages(b_idx, Lambda[run-1], Lambda[run-1][b_idx], m2f_theta_trimmed[l_idx])

        # Append processed data to the respective lists
        timestamp.append(B_t[run-1][b_idx])  # [s]
        wavelength.append(Lambd)  # [m]
        ellipticity.append(ep)  # [rad]
        angle.append(th)  # [rad]
        m2f_angle.append(m2f_th)  # [rad]
        frequency = y_weightedmean  # [Hz]
        frequency_sem = y_residualstd / np.sqrt(len(freq))  # [Hz]
        
        x = timestamp[0] / 60 # [min]
        y = m2f_angle[0] # [rad]

        # Fitting the data to a linear model
        y_fit, slope, intercept, m2f_weightedmean, residual, residualstd = self.analyzer.drift_fit(x, y, 3)
        y_linearfit, fitted_k, fitted_b, sigma_k, sigma_b, residuals, m2f_std, chi2_value, dof = self.analyzer.linear_fit(x, y, slope, intercept, y_residualstd)
        m2f_sem = m2f_std / np.sqrt(len(y))

        print("Average frequency detuning: {:.2f} MHz".format(frequency*1e-6 - self.consts.K39_D2_Hz * 1e-6))
        print("Standard error of the mean of frequency detuning: {:.2f} MHz".format(frequency_sem*1e-6))
        print("Average Faraday rotation: {:.2f} μrad".format(m2f_weightedmean*1e6))
        print("Standard error of the mean of Faraday rotation: {:.0f} nrad".format(m2f_sem*1e9))

        return frequency, frequency_sem, m2f_weightedmean, m2f_sem
    
    def process_B_field(self, gaussmeter_path, run):
        """
        Calculate the magnetic field strength and its variation.
        :param gaussmeter_path: Path to Gaussmeter data
        :param run: Current run index
        :return: Average magnetic field strength and its variation
        """
        # Read Gaussmeter data
        timestamps, B0s, temps = self.reader.read_gaussmeter(gaussmeter_path)

        # Instrument uncertainties and resolution
        instrument_uncertainty = 0.0005
        instrument_resolution = 0.02*1e-3

        # Compute the arithmatic mean and standard error of the mean of temperature measurements
        T_mean = np.mean(temps[run-1])
        T_sem = np.std(temps[run-1], ddof=1) / np.sqrt(len(temps[run-1]))

        # Linear fit of magnetic field measurements to estimate the drift and standard deviation
        B_linearfit, B_slope, B_intercept, B_weightedmean, B_residuals, B_residualstd = self.analyzer.drift_fit(timestamps[run-1], B0s[run-1], 20)

        # Calculate the uncertainties
        linear_sem = B_residualstd / np.sqrt(len(B0s[run-1]))
        sem1 = max(B_weightedmean * instrument_uncertainty, linear_sem)

        # Drift sine fit of magnetic field measurements to estimate the drift and standard deviation
        k0 = B_slope # slope guess
        b0 = B_intercept # intercept guess
        a0 = B_residualstd * 1e-3 # amplitude guess
        phi0 = 0 # phase guess  
        y_sinefit, fitted_k, fitted_b, B_amplitude, fitted_phi, sigma_k, sigma_b, sigma_a, sigma_phi, residuals, B_std, chi2_value, dof = self.analyzer.drift_sine_fit(timestamps[run-1], B0s[run-1], k0, b0, a0, phi0, sem1)
        B_sem = B_std / np.sqrt(len(B0s[run-1]))

        print("Average magnetic field strength: {:.3f} G".format(B_weightedmean))
        print("Magnetic field modulation amplitude: {:.1f} mG".format(B_amplitude*1e3))
        print("Standard error of the mean of magnetic field strength: {:.1f} mG".format(B_sem*1e3))
        return B_weightedmean, B_amplitude, B_sem, T_mean, T_sem
    
    def write_to_file(self, gaussmeter_path, lambda_path, lockin_path, path, filename, dtype, n, run):
        """
        Write the processed data to csv file
        """
        frequency, frequency_std, m2f_weightedmean, m2f_sem = self.process_physics(lambda_path, lockin_path, dtype, n, run)
        B_weightedmean, B_amplitude, B_sem, T_mean, T_sem = self.process_B_field(gaussmeter_path, run)
        data = [frequency, frequency_std, m2f_weightedmean, m2f_sem, -B_weightedmean, B_amplitude, B_sem, T_mean]

        try:
            file_path = os.path.join(path, filename)
            header = 'Frequency (Hz), Frequency_sem (Hz), Faraday_rotation (rad), Faraday_rotation_sem (rad), Longitudinal_magnetic_field (G), Magnetic_field_modualtion_amplitude (G), Magnetic_field_sem (G), Temperature (C)\n'
            if not os.path.isfile(file_path):
                with open(file_path, "w") as file:
                    file.write(header)

            # Check for duplicate data
            with open(file_path, "a+") as file:
                file.seek(0)
                existing_data = file.readlines()
                for line in existing_data:
                    existing_columns = line.strip().split(",")
                    new_columns = [str(item) for item in data]
                    if existing_columns == new_columns:
                        print("Duplicate data detected, abort writing.")
                        return
            
            # Append new data
            with open(file_path, "a") as file:
                file.write(",".join(map(str, data)) + "\n")

            print("Data appended to the file successfully.")
        except Exception as e:
            print(f"An error occurred while saving data to the file: {e}")

    def plot_triple_modulation_data(self, path):
        df = pd.read_csv(path, sep=',', header=None, skiprows=1,
                                 names=['Frequency (Hz)', 'Frequency_sem (Hz)', 
                                        'Faraday_rotation (rad)', 'Faraday_rotation_sem (rad)', 
                                        'Longitudinal_magnetic_field (G)', 'Magnetic_field_modualtion_amplitude (G)', 
                                        'Magnetic_field_sem (G)', 'Temperature (C)'])
    
        # Group data into chunks of 5 rows each
        group_size = 5
        num_groups = len(df) // group_size

        # Initialize lists to store results
        x_avg_list = []
        x_err_list = []
        y_avg_list = []
        y_err_list = []
        detun_list = []
        freq_err_list = []

        # Process each group
        for i in range(num_groups-1):
            group = df.iloc[i*group_size : (i+1)*group_size]

            # Frequency (x-axis) calculations
            freq_vals = group['Frequency (Hz)'].values
            freq_sem = group['Frequency_sem (Hz)'].values

            # Calculate frequency statistics
            f_unitfactor = 1e-6 # Convert Hz to MHz
            freq_avg = np.mean(freq_vals * f_unitfactor)
            var_between_freq = np.std(freq_vals * f_unitfactor, ddof=1)**2
            avg_within_var_freq = np.mean((freq_sem * f_unitfactor)**2)
            total_var_freq = (var_between_freq + avg_within_var_freq) / group_size
            f_err = np.sqrt(total_var_freq)

            # Magnetic field (x-axis) calculations
            b_vals = group['Magnetic_field_modualtion_amplitude (G)'].values
            b_sem = group['Magnetic_field_sem (G)'].values
            
            # Calculate x statistics
            b_unitfactor = 1e3
            x_avg = np.mean(b_vals * 2 * b_unitfactor)
            var_between_x = np.std(b_vals * 2 * b_unitfactor, ddof=1)**2
            avg_within_var_x = np.mean((b_sem * b_unitfactor)**2)
            total_var_x = (var_between_x + avg_within_var_x) / group_size
            x_err = np.sqrt(total_var_x)
            
            # Faraday rotation (y-axis) calculations
            fr_vals = group['Faraday_rotation (rad)'].values
            fr_sem = group['Faraday_rotation_sem (rad)'].values
            
            # Calculate y statistics
            fr_unitfactor = 1e6
            y_avg = np.mean(fr_vals * fr_unitfactor)
            var_between_y = np.std(fr_vals * fr_unitfactor, ddof=1)**2
            avg_within_var_y = np.mean((fr_sem * fr_unitfactor)**2)
            total_var_y = (var_between_y + avg_within_var_y) / group_size
            y_err = np.sqrt(total_var_y)
            
            # Store results
            x_avg_list.append(x_avg)
            x_err_list.append(x_err)
            y_avg_list.append(y_avg)
            y_err_list.append(y_err)
            detun_list.append(freq_avg - self.consts.K39_D2_Hz * 1e-6)
            freq_err_list.append(f_err)

        # Convert to numpy arrays
        x_data = np.array(x_avg_list)
        x_err = np.array(x_err_list)
        y_data = np.array(y_avg_list)
        y_err = np.array(y_err_list)
        detun_data = np.array(detun_list)
        freq_err = np.array(freq_err_list)

        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        y_fit, slope, intercept, y_weightedmean, residual, y_residualstd = self.analyzer.drift_fit(x_data, y_data, 3)
        y_linearfit, fitted_k, fitted_b, sigma_k, sigma_b, residuals, std, chi_value, dof = self.analyzer.linear_fit(x_data, y_data, slope, intercept, y_residualstd)

        # Shot-noise-limited angular sensitivity
        spectral_responsivity = 0.56 # A/W
        QE = spectral_responsivity * self.consts.h * self.consts.c / (self.consts.e * 766.7e-9) # Quantum efficiency
        snl_sensitivity = np.sqrt(self.consts.h * self.consts.c / (2 * QE * 766.7e-9 * 62.5e-6)) # [rad/sqrt(Hz)]
        print('Shot-noise-limited angular sensitivity = {:.2f} x 1e-8 rad/sqrt(Hz)'.format(snl_sensitivity*1e8))

        # Calculate the x-axis intercept
        x_intercept = -fitted_b / fitted_k
        x_fit = np.linspace(x_intercept, max(x_data), 6)
        y_fit = fitted_k * x_fit + fitted_b
        chi2 = np.sum((residuals / y_err)**2)  # Raw chi-squared
        print('Average ferquency detuning = {:.2f} MHz'.format(np.mean(detun_list)))
        print('Detuning error = {:.2f} MHz'.format(np.mean(freq_err)))
        print('Slope = {:.2f} μrad/mG'.format(fitted_k))
        print('Intercept = {:.2f} μrad'.format(fitted_b))
        print('Minimum detectable magnetic field variation:', x_intercept, 'mG')
        print('Chi-squared / dof:', chi2, '/', dof)

        y_lower_bound = y_fit - std
        y_upper_bound = y_fit + std

        ax.errorbar(x_data, y_data, xerr=x_err, yerr=y_err, fmt='.', capsize=5, color='b', label='Measured Faraday Rotation')
        ax.plot(x_fit, y_fit, 'r--', label=f'Linear Fit')
        ax.fill_between(x_fit, y_lower_bound, y_upper_bound, facecolor='C3', alpha=0.4, label=f'$\\sigma$={round(std,2):.2f} μrad')
        ax.set_xlabel(r'Longitudinal Magnetic Field Variation, $\Delta B_z$ (mG)', fontsize=30)
        ax.set_ylabel(r'Faraday Rotation Variation, $\Delta\theta$ (μrad)', fontsize=30)
        ax.tick_params(axis='both', which='major', labelsize=30)
        ax.tick_params(axis='both', which='minor', labelsize=30)
        # ax.set_title(f'$\chi^2/\\text{{dof}}$={chi2:.1f}/{dof}', fontsize=25)
        # ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=30)
        plt.tight_layout()
        save_path = os.path.join(plots, "Faraday_rotation_vs_magnetic_field_modulation(tightlayout).png")
        plt.savefig(save_path)
        plt.show()

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
    gaussmeter = os.path.join(K_vapor, 'Gaussmeter_data')
    lockins = os.path.join(K_vapor, 'Lockins_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data', 'Triple_modulation_measurements')

    faraday = FaradayRotation()
    date_input = '02-26-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    wavelengthmeter_path = glob.glob(os.path.join(wavelengthmeter, date, '*.csv'))
    gaussmeter_path = glob.glob(os.path.join(gaussmeter, date, '*.csv'))
    lockins_path = glob.glob(os.path.join(lockins, date, '*.lvm'))
    filename = 'Processed_tri_mod_data.csv'

    faraday.plot_triple_modulation_data(os.path.join(processed_path, filename))

    # for i in range(5, 10):
        # faraday.write_to_file(gaussmeter_path, wavelengthmeter_path, lockins_path, processed_path, filename, 'X', 10, i)