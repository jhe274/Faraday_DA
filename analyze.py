import numpy as np
from read import Read
import scipy.special
from scipy.optimize import curve_fit

class Analyze:

    def __init__(self):
        self.reader = Read()

    def FR_double_Kvapor(self, lockins_path):
        """
        Analyzed data from double modulated measurements
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockins_path)
        R1f, R2f, Rdc, epsilon, theta = [], [], [], [], []
        for i in range(len(lockins_path)):
            R1f.append(np.sqrt(X1f[i] ** 2 + Y1f[i] ** 2))                                                                  # 1f Magnitude: [V]
            R2f.append(np.sqrt(X2f[i] ** 2 + Y2f[i] ** 2))                                                                  # 2f Magnitude: [V]
            Rdc.append(np.sqrt(Xdc[i] ** 2 + Ydc[i] ** 2))                                                                  # dc Magnitude: [V]
            epsilon.append(R1f[i] / ( 2 * np.pi * scipy.special.jv(1,2.405) * Rdc[i]))                                      # Ellipticity: [rad]
            theta.append(R2f[i] / (2 * np.pi * scipy.special.jv(2,2.405) * Rdc[i] * np.sqrt(1 - 4 * epsilon[i]**2)))        # Rotation: [rad]
        
        return para, lockins_t, R1f, R2f, Rdc, epsilon, theta

    def FR_triple_Kvapor(self, lockins_path):
        """
        Analyzed data from triple modulated measurements
        """
        para, Lockins_t, Xmod, Ymod, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockins_path)
        Rmod, R1f, R2f, Rdc, epsilon, theta = [], [], [], [], [], []
        for i in range(len(lockins_path)):
            Rmod.append(np.sqrt(Xmod[i] ** 2 + Ymod[i] ** 2))                                                               # mod Magnitude: [V]
            R1f.append(np.sqrt(X1f[i] ** 2 + Y1f[i] ** 2))                                                                  # 1f Magnitude: [V]
            R2f.append(np.sqrt(X2f[i] ** 2 + Y2f[i] ** 2))                                                                  # 2f Magnitude: [V]
            Rdc.append(np.sqrt(Xdc[i] ** 2 + Ydc[i] ** 2))                                                                  # dc Magnitude: [V]
            epsilon.append(R1f[i] / (2 * np.pi * scipy.special.jv(1,2.405) * Rdc[i]))                                       # Ellipticity: [rad]
            theta.append(np.sqrt(2) * para[i][3] * Rmod[i] / 
                         (10 * np.pi * scipy.special.jv(2,2.405) * Rdc[i] * np.sqrt(1 - 4 * epsilon[i]**2)))                # Rotation: [rad]
        
        return para, Lockins_t, Rmod, R2f, Rdc, theta

    def check_calib(self, lambd, theta):
        '''
        Check the timestamp when Bristol starts and ends its self-calibration as calib
        and the corresponding timestamp in lock-ins as idx1, idx2
        '''
        # Bristol wavelength meter has several ways(manual, time, temperature) for self-calibration
        # every time the instrument calibrates there will be a ~2s blank in the wavelength measurements
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
        if bristol_t[-1] > lock_in_t[-1]:
            print(lambd[-1], lock_in_t[-1])
            end_t = np.argmin(np.abs(bristol_t - lock_in_t[-1]))
            trimmed_bristol_t, trimmed_lambd, trimmed_lock_in_t, trimmed_theta = bristol_t[:end_t+1], lambd[:end_t+1], lock_in_t, theta
        else:
            end_t = np.argmin(np.abs(lock_in_t - bristol_t[-1]))
            # print(lock_in_t)
            # print(theta)
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
