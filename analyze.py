import numpy as np
from read import Read
import scipy.special
from scipy.optimize import curve_fit

class Analyze:

    def __init__(self):
        self.reader = Read()

    def R_lockins(self, lockins_path):
        """
        Calculate R values from lock-in amplifiers
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockins_path)
        R1f, R2f, Rdc = [], [], []
        for i in range(len(lockins_path)):
            # Rmod.append(np.sqrt(Xmod[i] ** 2 + Ymod[i] ** 2))                                                               # mod Magnitude: [V]   
            R1f.append(np.sqrt(X1f[i] ** 2 + Y1f[i] ** 2))                                                                  # 1f Magnitude: [V]
            R2f.append(np.sqrt(X2f[i] ** 2 + Y2f[i] ** 2))                                                                  # 2f Magnitude: [V]
            Rdc.append(np.sqrt(Xdc[i] ** 2 + Ydc[i] ** 2))                                                                  # dc Magnitude: [V]

        return para, lockins_t, R1f, R2f, Rdc
        # return para, lockins_t, Rmod, R1f, R2f, Rdc

    def FR_double_Kvapor(self, lockins_path, S1f, S2f, Sdc):
        """
        Analyzed data from double modulated measurements
        """
        epsilon, theta = [], []
        for i in range(len(lockins_path)):
            epsilon.append(S1f[i] / ( 2 * np.pi * scipy.special.jv(1,2.405) * Sdc[i]))                                      # Ellipticity: [rad]
            theta.append(S2f[i] / (2 * np.pi * scipy.special.jv(2,2.405) * Sdc[i] * np.sqrt(1 - 4 * epsilon[i]**2)))        # Rotation: [rad]
        
        return epsilon, theta

    def FR_triple_Kvapor(self, lockins_path, Smod, S1f, S2f, Sdc):
        """
        Analyzed data from triple modulated measurements
        """
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self. reader.lockins(lockins_path)
        # para, lockins_t, Rmod, R1f, R2f, Rdc = self.R_lockins(lockins_path)
        epsilon, theta = [], []
        for i in range(len(lockins_path)):
            epsilon.append(S1f[i] / (2 * np.pi * scipy.special.jv(1,2.405) * Sdc[i]))                                       # Ellipticity: [rad]
            theta.append(np.sqrt(2) * para[i][3] * Smod[i] / 
                         (10 * np.pi * scipy.special.jv(2,2.405) * Sdc[i] * np.sqrt(1 - 4 * epsilon[i]**2)))                # Rotation: [rad]
        
        return epsilon, theta

    def check_calib(self, lambd, theta):
        '''
        Check the timestamp when Bristol starts and ends its self-calibration as calib
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
    
    def filter_data(self, t, lambd):
        '''
        Filtering lambda values based on bounds from Bristol measurements
        '''
        try:
            lambda_ubound = 766.72 * 1e-9
            lambda_lbound = 766.68 * 1e-9
            condition = np.where(lambd < lambda_ubound)
            # condition = np.logical_and(lambd > lambda_lbound, lambd < lambda_ubound)
            filtered_t = t[condition]
            filtered_lambd = lambd[condition]
            if filtered_lambd.size == 0:
                raise ValueError("No lambda values found within the specified bounds.")
        except ValueError as e:
            print(e)
            t, lambd = t, lambd

        return filtered_t, filtered_lambd
    
    def trim_data(self, x1, y1, x2, y2):
        '''
        Trimming data based on the shorter measurement from either Bristol or lock-ins
        '''
        if x1[-1] > x2[-1]:
            print(y1[-1], x2[-1])
            idx = np.argmin(np.abs(x1 - x2[-1]))
            trim_x1, trim_y1, trim_x2, trimmed_y2 = x1[:idx+1], y1[:idx+1], x2, y2
        else:
            idx = np.argmin(np.abs(x2 - x1[-1]))
            trim_x1, trim_y1, trim_x2, trim_y2 = x1, y1, x2[:idx+1], y2[:idx+1]

        return trim_x1, trim_y1, trim_x2, trim_y2

    def calculate_interval_and_indices(self, x1, x2, TC, n):
        '''
        Timestamp matching between lock-ins and Bristol
        '''
        # Data acquisition using ReadFaradayLockinsRealtimeplot.vi typically
        # has a update rate about 10Hz (100ms) per data point, and similarly using lock-ins(x2) internal buffer
        # the storage interval settings are typically slower than Bristol(x1), where Bristol has
        # a minimum frame rate 20Hz (50ms) per data point
        if TC > 200e-3:
            interval = np.arange(TC * n, x2[-1] + TC, TC)
            x2_idx = np.searchsorted(x2, interval, side='left')[:-1]
            x1_idx = np.searchsorted(x1, x2[x2_idx], side='left')
        else:
            x2_idx = np.arange(n,len(x2)-1,1)
            x1_idx = np.searchsorted(x1, x2[x2_idx], side='left')

        # Removing duplicate timestamps from x1_idx, and their corresponding timestamps in x2
        uni_idx = np.unique(x1_idx, return_index=True)[1]
        dup_idx = np.setdiff1d(np.arange(len(x1_idx)), uni_idx)
        l_idx = np.delete(x2_idx, dup_idx)
        b_idx = np.unique(x1_idx)

        return l_idx, b_idx
    
    def calculate_averages(self, b_idx, y1, filtered_y1, y2):
        '''
        Average function to match lock-ins readings and wavelength measurements
        '''
        calib, theta = self.check_calib(filtered_y1, y2)
        averages = []
        for i in np.arange(len(b_idx) - 1):
            if calib.size > 0 and i + 1 >= calib[0] and i + 1 < calib[-1] + 2:
                pass
            else:
                segment = y1[b_idx[i]:b_idx[i + 1]]
                avg_value = np.mean(segment)
                averages.append(avg_value)

        return np.array(averages), theta
    
    def convert_to_float_array(self, arr):
        '''
        Convert processed data from object to float array
        '''
        return np.array([np.array(lst, dtype=np.float64) for lst in np.array(arr, dtype=object)], dtype=object)
    
    def smooth(self, data, width=100):
        len_data = len(data)
        smoothed_data = np.zeros(len_data)
        for i in range(len_data):
            lower = max(0,i-width//2)
            upper = min(i+width//2, len_data)
            avg = np.mean(data[lower:upper])
            smoothed_data[i] = avg

        return smoothed_data
