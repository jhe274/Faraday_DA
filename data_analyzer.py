import numpy as np
from data_reader import DataReader as Read
import scipy.special
from scipy.optimize import curve_fit

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
            tuple: (para, lockins_t, R1f, R2f, Rdc, Rmod) with R values in volts.
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xmod, Ymod = self.reader.read_lockins(lockins_path)

        R1f = [np.hypot(X1f[i], Y1f[i]) for i in range(len(lockins_path))]
        R2f = [np.hypot(X2f[i], Y2f[i]) for i in range(len(lockins_path))]
        Rdc = [np.hypot(Xdc[i], Ydc[i]) for i in range(len(lockins_path))]
        Rmod = [np.hypot(Xmod[i], Ymod[i]) for i in range(len(lockins_path))]

        return para, lockins_t, R1f, R2f, Rdc, Rmod
    
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
    
    def modulated_angle(self, lockins_path, V1f, Vdc, Vmod):
        """
        Computes the modulated rotation angle.

        Parameters:
            lockins_path (list): List of lock-in amplifier file paths.
            V1f (list): First harmonic voltage component.
            Vdc (list): DC voltage component.
            Vmod (list): Modulation voltage component.

        Returns:
            tuple: (theta) in radians.
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xmod, Ymod = self.reader.read_lockins(lockins_path)
        epsilon, epsilon_approx = self.ellipticity(lockins_path, V1f, Vdc)
        print([para[i][3] for i in range(len(lockins_path))])
        theta_mod = [
            0.5 * np.arcsin(
                np.clip(np.sqrt(2) * para[i][3] * Vmod[i] / (
                5 * np.pi * scipy.special.jv(2, 2.405) * Vdc[i] * np.sqrt(1 - 4 * epsilon_approx[i]**2)), -1, 1)
            ) 
            for i in range(len(lockins_path))
        ]

        return theta_mod
    
    def absorbance_difference(self, epsilon, l):
        """Computes absorbance difference."""
        return 4 * np.array(epsilon) / l

    def refractive_indices_difference(self, theta, l, wavelength):
        """Computes refractive indices difference."""
        return np.array(wavelength) * np.array(theta) / (np.pi * l)

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
            tuple: (filtered_t, filtered_x)
        """
        x_ubound, x_lbound = 766.72e-9, 766.68e-9
        condition = np.logical_and(x > x_lbound, x < x_ubound)

        filtered_t = t[condition]
        filtered_x = x[condition]

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

    def smooth(self, data, width=100):
        """Applies a simple moving average smoothing filter."""
        len_data = len(data)
        return np.array([np.mean(data[max(0, i-width//2):min(i+width//2, len_data)]) for i in range(len_data)])
