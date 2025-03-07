import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from data_reader import DataReader as Read
import scipy.special
from scipy.stats import linregress, chisquare
from scipy.signal import savgol_filter, find_peaks, butter, filtfilt
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from lmfit import Model

class DataAnalyzer:

    def __init__(self):
        """Initializes the Analyze class with a DataReader instance."""
        self.reader = Read()

    def R_lockins(self, lockins_path):
        """
        Computes R values from lock-in amplifier measurements.
        
        Parameters:
            lockins_path (list): List of lock-in amplifier file paths.
        
        Returns:
            tuple: (para, lockins_t, R1f, R2f, Rdc, Rm2f) with R values in volts.
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f = self.reader.read_lockins(lockins_path)

        R1f = [np.hypot(X1f[i], Y1f[i]) for i in range(len(lockins_path))]
        R2f = [np.hypot(X2f[i], Y2f[i]) for i in range(len(lockins_path))]
        Rdc = [np.hypot(Xdc[i], Ydc[i]) for i in range(len(lockins_path))]
        Rm2f = [np.hypot(Xm2f[i], Ym2f[i]) for i in range(len(lockins_path))]

        return para, lockins_t, R1f, R2f, Rdc, Rm2f
    
    def ellipticity(self, lockins_path, V1f, Vdc):
        """
        Computes ellipticity and its approximation.
        
        Parameters:
            lockins_path (list): List of lock-in amplifier file paths.
            V1f (list): First harmonic voltage component.
            Vdc (list): DC voltage component.
        
        Returns:
            tuple: (epsilon, epsilon_approx) in radians.
        """
        epsilon = [0.5 * np.arctanh(np.clip(2 * V1f[i] / (np.pi * scipy.special.jv(1, 2.405) * Vdc[i]), -1, 1)) for i in range(len(lockins_path))]
        epsilon_approx = [V1f[i] / (np.pi * scipy.special.jv(1, 2.405) * Vdc[i]) for i in range(len(lockins_path))]
        
        return epsilon, epsilon_approx

    def angle(self, lockins_path, V1f, V2f, Vdc):
        """
        Computes polarization rotation angle.

        Parameters:
            lockins_path (list): List of lock-in amplifier file paths.
            V1f (list): First harmonic voltage component.
            V2f (list): Second harmonic voltage component.
            Vdc (list): DC voltage component.
        
        Returns:
            tuple: (theta) in radians.
        """
        epsilon, epsilon_approx = self.ellipticity(lockins_path, V1f, Vdc)
        theta = [
            0.5 * np.arcsin(
                np.clip(2 * V2f[i] / (np.pi * scipy.special.jv(2, 2.405) * Vdc[i] * np.sqrt(1 - 4 * epsilon_approx[i]**2)), -1, 1)
            )
            for i in range(len(lockins_path))
        ]

        return theta
    
    def modulated_angle(self, lockins_path, V1f, Vdc, Vm2f):
        """
        Computes the modulated rotation angle.

        Parameters:
            lockins_path (list): List of lock-in amplifier file paths.
            V1f (list): First harmonic voltage component.
            Vdc (list): DC voltage component.
            Vm2f (list): B-field modulation voltage component.

        Returns:
            tuple: (theta) in radians.
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f = self.reader.read_lockins(lockins_path)
        epsilon, epsilon_approx = self.ellipticity(lockins_path, V1f, Vdc)
        theta_m2f = [
            0.5 * np.arcsin(
                np.clip(np.sqrt(2) * para[i][3] * Vm2f[i] / (
                5 * np.pi * scipy.special.jv(2, 2.405) * Vdc[i] * np.sqrt(1 - 4 * epsilon_approx[i]**2)), -1, 1)
            ) 
            for i in range(len(lockins_path))
        ]

        return theta_m2f
    
    def absorbance_difference(self, epsilon, l):
        """Computes absorbance difference."""
        return 4 * np.array(epsilon) / l # [rad/m]

    def refractive_indices_difference(self, theta, l, wavelength):
        """Computes refractive indices difference."""
        return np.array(wavelength) * np.array(theta) / (l * np.pi) # [rad]

    def check_calib(self, x, y):
        """
        Detects and handles Bristol self-calibration.
        """
        calib = np.where(x == 0)[0]
        if calib.size > 0:
            print('Self-calibration detected in Bristol measurements...')
            y[calib[0]:calib[-1] + 2] = 0
        else:
            print('No self-calibration detected...')
        
        return calib, y
    
    def filter_data(self, t, x):
        """
        Filters data values within a specified range.
        
        Parameters:
            t (array-like): Time array.
            x (array-like): Data array (e.g., wavelength).
        
        Returns:
            np.ndarray, np.ndarray: Filtered time and data arrays.
        """
        t = np.asarray(t)
        x = np.asarray(x)

        # x_ubound, x_lbound = 766.701e-9, 766.6997e-9
        # condition = np.logical_and(x > x_lbound, x < x_ubound)

        # filtered_t = t[condition]
        # filtered_x = x[condition]

        filtered_t = t
        filtered_x = x

        if filtered_x.size == 0:
            raise ValueError("No values found within the specified bounds.")
        

        return filtered_t, filtered_x
    
    def trim_data(self, x1, y1, x2, y2):
        """Trims data to match shorter timestamps between two measurement sources."""
        if x1[-1] > x2[-1]: 
            idx = np.argmin(np.abs(x1 - x2[-1]))
            return x1[:idx+1], y1[:idx+1], x2, y2
        else:  
            idx = np.argmin(np.abs(x2 - x1[-1]))
            return x1, y1, x2[:idx+1], y2[:idx+1]

    def calculate_interval_and_indices(self, x1, x2, TC, n):
        """Matches timestamps between lock-ins and Bristol data."""
        if TC > 200e-3:
            interval = np.arange(TC * n, x2[-1] + TC, TC)
            x2_idx = np.searchsorted(x2, interval, side='left')[:-1]
        else:
            x2_idx = np.arange(n, len(x2) - 1, 1)

        x1_idx = np.searchsorted(x1, x2[x2_idx], side='left')

        # Removing duplicate timestamps from x1_idx, and their corresponding timestamps in x2
        uni_idx = np.unique(x1_idx, return_index=True)[1]
        dup_idx = np.setdiff1d(np.arange(len(x1_idx)), uni_idx)
        l_idx = np.delete(x2_idx, dup_idx)
        b_idx = np.unique(x1_idx)

        return l_idx, b_idx
    
    def calculate_averages(self, indices, x, filtered_x, y):
        """Computes moving averages for indexed data while handling self-calibration events."""
        calib, filtered_y = self.check_calib(filtered_x, y)
        averages = [
            np.mean(x[indices[i]:indices[i + 1]]) 
            for i in np.arange(len(indices) - 1) 
            if calib.size == 0 or not (calib[0] <= i + 1 < calib[-1] + 2)
        ]
        return np.array(averages), filtered_y
    
    def bin_data(self, data, bin_size):
        """
        Bins data by averaging over specified bin_size.
        
        Parameters:
            data (array-like): Input data array.
            bin_size (int): Number of data points per bin.
        
        Returns:
            np.ndarray: Binned data array.
        """
        num_bins = len(data) // bin_size
        binned_data = np.array([np.mean(data[i * bin_size:(i + 1) * bin_size]) for i in range(num_bins)])
        
        # Handle remaining points
        if len(data) % bin_size != 0:
            remaining_mean = np.mean(data[num_bins * bin_size:])
            binned_data = np.append(binned_data, remaining_mean)

        return binned_data
    
    def rolling_variance(self, data, window):
        """
        Compute the rolling variance of a 1D array, ignoring NaNs.

        Parameters:
        data (numpy.ndarray): Input data array.
        window (int): The size of the rolling window.

        Returns:
        numpy.ndarray: Array of rolling variances.
        """
        if window < 1:
            raise ValueError("Window size must be at least 1.")
        if window > len(data):
            raise ValueError("Window size must not be larger than the data length.")

        # Handle NaNs in input data (replace with local mean if needed)
        data = np.where(np.isnan(data), np.nanmean(data), data)

        # Use sliding_window_view to create rolling windows
        rolling_windows = sliding_window_view(data, window)

        # Compute variance for each window
        variances = np.nanvar(rolling_windows, axis=-1, ddof=1)  # ddof=1 for sample variance

        # Pad the result with the median variance instead of NaN
        pad_width = window - 1
        median_variance = np.nanmedian(variances)  # Compute median to replace NaNs
        variances = np.pad(variances, (pad_width, 0), mode='constant', constant_values=median_variance)

        return variances

    def moving_average(self, y, bin):
        """
        ✅ Pros: Simple and effective for removing high-frequency noise.
        ❌ Cons: Can distort peak shapes and shift the spectrum.
        """
        return np.convolve(y, np.ones(bin)/bin, mode='same')
    
    def polynomial_smoothing(self, y, length, order):
        """"
        ✅ Pros: Preserves peak shapes and spectral features.
        ❌ Cons: Not ideal for extremely noisy data or very small datasets.
        """
        return savgol_filter(y, window_length=length, polyorder=order)
    
    def gaussian_smoothing(self, y, sigma):
        """
        ✅ Pros: Good balance between noise reduction and feature preservation.
        ❌ Cons: Can introduce small artifacts at spectrum edges.
        """
        return gaussian_filter1d(y, sigma)
    
    def bandpass_filter(self, data, lowcut, highcut, fs, order=4):
        nyquist = 0.5 * fs
        low = lowcut / nyquist
        high = highcut / nyquist
        b, a = butter(order, [low, high], btype='band', analog=False)

        return filtfilt(b, a, data)
    
    def fft_peak(self, x, y):
        N = len(y)
        dt = np.mean(np.diff(x))
        y_fft = np.fft.fft(y)
        y_fft_freq = np.fft.fftfreq(N, dt)

         # Take only positive frequencies
        half_N = N // 2  # Half-point index
        positive_freqs = y_fft_freq[:half_N]
        fft_magnitudes = np.abs(y_fft[:half_N])

        # Correct normalization for amplitude spectrum
        amplitude_spectrum = fft_magnitudes / N  # Normalization
        amplitude_spectrum[1:] *= 2  # Double non-DC components

        # Find peaks in the FFT spectrum
        peaks, properties = find_peaks(amplitude_spectrum, height=0.001)  # Adjust threshold if needed

        if len(peaks) > 0:
            # Find the closest peak to 0.5 Hz
            closest_peak_index = np.argmin(np.abs(positive_freqs[peaks] - 0.5))
            dominant_freq = positive_freqs[peaks][closest_peak_index]
            print(f"Detected Dominant Frequency Near 0.5 Hz: {dominant_freq:.6f} Hz")

            # Find the closest peak to 0.00086 Hz
            # closest_peak_index = np.argmin(np.abs(positive_freqs[peaks] - 0.00086))
            # dominant_freq = positive_freqs[peaks][closest_peak_index]
            # print(f"Detected Dominant Frequency Near 0.00086 Hz: {dominant_freq:.6f} Hz")
        else:
            dominant_freq = None
            print("No peaks detected in the specified range.")
        
        return positive_freqs, fft_magnitudes, amplitude_spectrum, dominant_freq
    
    def noise_spectral_density(self, x, y):
        N = len(y)
        dt = np.mean(np.diff(x))
        y_fft = np.fft.fft(y)
        y_fft_freq = np.fft.fftfreq(N, dt)

        positive_freqs = y_fft_freq[y_fft_freq >= 0]
        fft_magnitudes = np.abs(y_fft)[y_fft_freq >= 0]

        # Apply single-sided scaling (×2 for non-DC components)
        non_dc_mask = positive_freqs > 0     # Mask to exclude DC (0 Hz)
        nsd = fft_magnitudes.copy()
        nsd[non_dc_mask] *= 2                # Double non-DC amplitudes
        nsd /= np.sqrt(N / dt)               # Normalize to 1/√Hz

        return positive_freqs, nsd

    def noise_floor(self, x, y):
        y_fft_freq, y_fft, dominant_freq, amplitude_spectrum = self.fft_peak(x, y)
        # Compute power spectral density
        psd = amplitude_spectrum**2

        # Define noise floor as the median power at non-dominant frequencies
        noise_floor = np.median(psd[y_fft_freq > 0.1])

        # Identify signal-to-noise ratio (SNR)
        peak_power = psd[np.argmin(np.abs(y_fft_freq - dominant_freq))]
        snr = peak_power / noise_floor

        print(f"Estimated Noise Floor (dB): {10 * np.log10(noise_floor):.2f}")
        print(f"Signal Power (dB) at {dominant_freq:.6f} Hz: {10 * np.log10(peak_power):.2f}")
        print(f"Signal-to-Noise Ratio (dB): {10 * np.log10(snr):.2f}")

        return noise_floor, psd, peak_power
    
    def bootstrap_errors(self, x, y, func, initial_guess, num_samples=100):
        """Bootstrap resampling to estimate parameter uncertainties."""
        fitted_params = []
        for _ in range(num_samples):
            resampled_y = y + np.random.normal(0, np.std(y - func(x, *initial_guess)), size=len(y))
            try:
                popt, _ = curve_fit(func, x, resampled_y, p0=initial_guess)
                fitted_params.append(popt)
            except RuntimeError:
                continue  # Skip failed fits

        fitted_params = np.array(fitted_params)
        param_std = np.std(fitted_params, axis=0)  # Standard deviation as error estimate
        return param_std

    def drift_sine_wave(self, t, k, b, a, phi):
        freq = 0.5
        # freq = 0.000867 * 60
        return b + k * t + a * np.sin(2 * np.pi * freq * t + phi)

    def line(self, x, a, b):
        return a * x + b
    
    def drift_fit(self, x, y, window):
        # Perform linear regression
        slope, intercept, _, _, _ = linregress(x, y)
        
        # Compute fitted values
        y_fit = slope * x + intercept

        # Compute the residuals
        residuals = y - y_fit

        # Compute the weights
        residuals_variance = self.rolling_variance(residuals, window)
        weights = 1 / residuals_variance

        # Compute the weighted mean
        y_weightedmean = np.average(y, weights=weights)

        # Compute standard deviation of residuals
        y_residualstd = np.std(residuals, ddof=1)

        return y_fit, slope, intercept, y_weightedmean, residuals, y_residualstd

    def linear_fit(self, x, y, k, b, sigma):
        """Curve fit using Scipy.optimize.curve_fit"""
        # initial guess
        initial_guess = [k, b]
        # Fit a sine wave to the data
        popt, pcov = curve_fit(self.line, x, y, p0=initial_guess)
        # Extract the fitted parameters
        fitted_k, fitted_b = popt
        sigma_k = np.sqrt(pcov[0, 0])
        sigma_b = np.sqrt(pcov[1, 1])
        # Print the true and fitted parameters
        print(f"Fitted Parameters:")
        print(f"Slope = {fitted_k}")
        print(f"Intercept = {fitted_b}")
        # Generate the fitted sine wave
        y_linearfit = self.line(x, *popt)
        dof = len(y) - len(popt)
        residuals = y - y_linearfit
        residual_std = np.std(residuals, ddof=len(popt))
        # sigma = self.bootstrap_errors(x, y, self.line, initial_guess)
        chi2_value = np.sum((residuals / sigma) ** 2)
        print(f"Reduced Chi-squared: {chi2_value:.1f}/{dof}")

        return y_linearfit, fitted_k, fitted_b, sigma_k, sigma_b, residuals, residual_std, chi2_value, dof
    
    def drift_sine_fit(self, x, y, k, b, a, phi, sigma):
        """Curve fit using Scipy.optimize.curve_fit"""
        # initial guess
        initial_guess = [k, b, a, phi]
        # Fit a sine wave to the data
        popt, pcov = curve_fit(self.drift_sine_wave, x, y, p0=initial_guess)
        # Extract the fitted parameters
        fitted_k, fitted_b, fitted_a, fitted_phi = popt
        sigma_k = np.sqrt(pcov[0, 0])
        sigma_b = np.sqrt(pcov[1, 1])
        sigma_a = np.sqrt(pcov[2, 2])
        sigma_phi = np.sqrt(pcov[3, 3])
        # Print the true and fitted parameters
        print(f"Fitted Parameters:")
        print(f"Slope = {fitted_k}")
        print(f"Intercept = {fitted_b}")
        print(f"Amplitude = {fitted_a}")
        print(f"Phase = {fitted_phi}")
        # Generate the fitted sine wave
        y_sinefit = self.drift_sine_wave(x, *popt)
        dof = len(y) - len(popt)
        residuals = y - y_sinefit
        residual_std = np.std(residuals, ddof=len(popt))
        # sigma = self.bootstrap_errors(x, y, self.drift_sine_wave, initial_guess)
        chi2_value = np.sum((residuals / sigma) ** 2)
        print(f"Reduced Chi-squared: {chi2_value:.1f}/{dof}")

        return y_sinefit, fitted_k, fitted_b, fitted_a, fitted_phi, sigma_k, sigma_b, sigma_a, sigma_phi, residuals, residual_std, chi2_value, dof