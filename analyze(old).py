import math
import numpy as np
from read import Read
import scipy.special
from scipy.optimize import curve_fit

class Analyze:

    def __init__(self):
        self.Read = Read()

    def Double_mod_theta(self, lock_in_path):
        # Faraday rotation angles of double modulated measurements
        para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod = self.Read.Lock_ins(lock_in_path)
        Rdc, R2f, theta = [], [], []

        for i in range(len(lock_in_path)):
            Rdc.append(np.sqrt(Xdc[i] ** 2 + Ydc[i] ** 2)) # [V]
            R2f.append(np.sqrt(X2f[i] ** 2 + Y2f[i] ** 2)) # [V]
            theta.append(R2f[i] / (math.pi * scipy.special.jv(2,2.405) * Rdc[i])) # [rad]
        
        return para, Lock_in_t, theta

    
    def Triple_mod_theta(self, lock_in_path):
        # Faraday rotation angles of triple modulated measurements
        para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod = self.Read.Lock_ins(lock_in_path)

        Rdc, R2f, Rmod, theta = [], [], [], []
        for i in range(len(lock_in_path)):
            Rdc.append(np.sqrt(Xdc[i] ** 2 + Ydc[i] ** 2)) # [V]
            R2f.append(np.sqrt(X2f[i] ** 2 + Y2f[i] ** 2)) # [V]
            Rmod.append(np.sqrt(Xmod[i] ** 2 + Ymod[i] ** 2)) # [V]
            theta.append(np.sqrt(2) * para[i][3] * Rmod[i]/(10 * math.pi * scipy.special.jv(2,2.405) * Rdc[i])) # [rad]
        
        return para, Lock_in_t, theta
    
    
    def check_calib(self, lambd, theta):
        '''
            Check the timestamp when Bristol starts and ends its self-calibration as calib
            and the corresponding timestamp in lock-ins as idx1, idx2
        '''
        # Bristol wavelength meter has several ways(manual, time, temperature) for self-calibration, every time the instrument calibrates there will be a ~2s blank in the wavelength measurements
        calib = np.where(lambd == 0)[0]
        if calib.size > 0:
            print('Self-calibration in Bristol measurements detected...')
            theta[calib[0]:calib[-1] + 2] = 0
        else:
            print('No self-calibration in Bristol measurements detected...')
            pass

        return calib, theta
    
    def filter_data(self, bristol_t, lambd):
        '''
            Filtering lambda values based on bounds
        '''
        try:
            lambda_ubound = 766.72 * 1e-9
            lambda_lbound = 766.68 * 1e-9

            condition = np.where(lambd < lambda_ubound)
            # condition = np.logical_and(lambd > lambda_lbound, lambd < lambda_ubound)
            filtered_bristol_t = bristol_t[condition]
            filtered_lambd = lambd[condition]

            if filtered_lambd.size == 0:
                raise ValueError("No lambda values found within the specified bounds.")
            
        except ValueError as e:
            print(e)
            bristol_t, lambd = bristol_t, lambd

        return filtered_bristol_t, filtered_lambd
    
    def trim_data(self, bristol_t, lambd, lock_in_t, theta):
        '''
            Trimming lock-in or Bristol data based on the shorter measurement
        '''
        if lambd[-1] > lock_in_t[-1]:
            end_t = np.argmin(np.abs(bristol_t - lock_in_t[-1]))
            trimmed_bristol_t, trimmed_lambd, trimmed_lock_in_t, trimmed_theta = bristol_t[:end_t+1], lambd[:end_t+1], lock_in_t, theta

        else:
            end_t = np.argmin(np.abs(lock_in_t - bristol_t[-1]))
            trimmed_bristol_t, trimmed_lambd, trimmed_lock_in_t, trimmed_theta = bristol_t, lambd, lock_in_t[:end_t+1], theta[:end_t+1]

        return trimmed_bristol_t, trimmed_lambd, trimmed_lock_in_t, trimmed_theta

    def calculate_interval_and_indices(self, bristol_t, lock_in_t, time_const, n):
        '''
            Timestamp matching between lock-ins and Bristol
        '''
        # Data acquisition using ReadFaradayLockinsRealtimeplot.vi typically has a update rate about 100ms per data point, 
        # making it significantly slower than Bristol, where Bristol has a minimum frame rate 20Hz that is 50ms per data point
        if time_const > 200e-3:
            interval = np.arange(time_const * n, lock_in_t[-1] + time_const, time_const)
            lock_in_idx = np.searchsorted(lock_in_t, interval, side='left')[:-1]
            Bristol_idx = np.searchsorted(bristol_t, lock_in_t[lock_in_idx], side='left')
        else:
            lock_in_idx = np.arange(n,len(lock_in_t)-1,1)
            Bristol_idx = np.searchsorted(bristol_t, lock_in_t[lock_in_idx], side='left')

        # Removing duplicate timestamps from Bristol_idx, and their corresponding timestamps in lock-ins
        uni_idx = np.unique(Bristol_idx, return_index=True)[1]
        dup_idx = np.setdiff1d(np.arange(len(Bristol_idx)), uni_idx)
        l_idx = np.delete(lock_in_idx, dup_idx)
        b_idx = np.unique(Bristol_idx)

        return l_idx, b_idx
    
    def calculate_averages(self, Bristol_idx, lambd, filtered_lambd, theta):
        '''
            Average function to match lock-in readings and wavelength measurements
        '''
        calib, theta = self.check_calib(filtered_lambd, theta)
        averages = []

        for i in np.arange(len(Bristol_idx) - 1):
            if calib.size > 0 and i + 1 >= calib[0] and i + 1 < calib[-1] + 2:
                pass
            else:
                segment = lambd[Bristol_idx[i]:Bristol_idx[i + 1]]
                avg_value = np.mean(segment)
                averages.append(avg_value)

        return np.array(averages), theta
    
    def interp_theta(self, bristol_t, lock_in_t, theta):
        # Interpolate polarization rotation measurements to wavelength measurements
        itp_theta = np.interp(bristol_t, lock_in_t, theta)
        # Handle NaN and Inf values
        itp_theta = np.nan_to_num(itp_theta, nan=0.0, posinf=0)
        
        return itp_theta
    
    def find_peaks_idx(self, x):
        # Find the peaks of wavelength scan
        dx = np.diff(x)
        dx = np.where(dx > 0, 1, -1)
        ddx = np.diff(dx)
        i_peaks = (np.where(ddx != 0)[0] + 1,)
        
        return i_peaks
    
    def modified_sine_wave(self, t, A, f, phi, offset, drift_slope):
        # Define the modified sine wave function
        return A * np.sin(2 * np.pi * f * t + phi) + offset + drift_slope * t
    
    def alan_variance(self, data):
        # Calculation of Alan variance
        n = len(data)
        mean_diff = np.diff(data).mean()
        alan_var = np.sum((np.diff(data) - mean_diff) ** 2) / (2 * n)
        
        return alan_var
    
    def filtered_theta_and_lambda(self, lambda_path, lock_in_path, i, n):
        # Calculate filterd polarization rotation and wavelength measurements
        Bristol_t, Lambda = self.Read.read_Bristol(lambda_path)
        para, Lock_in_t, theta = self.Triple_mod_theta(lock_in_path)
        Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i] = self.filter_and_trim_data(Bristol_t[i], Lambda[i], Lock_in_t[i], theta[i])
        # itp_theta = self.interp_theta(Bristol_t[i], Lock_in_t[i], theta[i])
        Bristol_idx, lock_in_idx = self.calculate_interval_and_indices(Bristol_t[i], Lock_in_t[i], para[i][4], n)

        N = len(lock_in_idx)
        lambd = Lambda[i][Bristol_idx[0]:] # [m]
        thet = theta[i][lock_in_idx] # [rad]
        Lambda_ave = np.average(lambd) # [m]
        Lambda_std = np.std(lambd, ddof=1.5) # [m], unbiased sample STD
        Lambda_ste = Lambda_std / np.sqrt(N) # [m]
        Theta_ave = np.average(thet) # [rad]
        Theta_std = np.std(thet, ddof=1.5) # [rad], unbiased sample STD
        Theta_ste = Theta_std / np.sqrt(N) # [rad]

        params, covariance = curve_fit(self.modified_sine_wave, Bristol_t[i][Bristol_idx[0]:], lambd, p0=[Lambda_std*2, 10, 0, Lambda_ave, 0])
        A_fit, f_fit, phi_fit, offset_fit, drift_slope_fit = params
        print("Fitted Delta lambda is " + str(round(np.abs(A_fit)*2*1e9,6)) + ' [nm]')
        
        return N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste