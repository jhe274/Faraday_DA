import re, os, glob
import datetime as dt
import numpy as np
import pandas as pd

class Read:

    def sort_key(self, path):
        '''
            Extract the numeric part after the last underscore in the file name
        '''
        match = re.search(r'(\d+)(?=\.[^.]*$|$)', path)
        if match:
            return int(match.group(0))
        else:
            return path
        
    
    def correct_negative_time_stamp(self, t):
        '''
            Correct the negative time stamps in the lock-in data file
        '''
        if np.any(t < 0):
            t_idx = np.where(t < 0)[0][0]
            t[t_idx:] = 24 * 60 * 60 + t[t_idx:]

    
    def Bristol(self, path):
        '''
            read Wavelength measurements from Brilstol 871 wavelength meter
        '''
        B_t, B_lambda = [], []

        for file in sorted(path, key=self.sort_key):
            df = pd.read_csv(file, sep=',', header=None, skiprows=1, 
                             names=['Timestamp', 'Instrument Status', 'Instrument Wavelength', 'Instrument Intensity'])
            
            df['Timestamp'] = pd.to_datetime(df['Timestamp'])
            start_time = df['Timestamp'].iloc[0]
            df['Timestamp'] = (df['Timestamp'] - start_time).dt.total_seconds()
            
            # Filter out error readings & mode hop readings
            # filtered_data = df[(df['Instrument Wavelength'] != 0) & (df['Instrument Wavelength'] < 767)]
    
            B_t.append(df['Timestamp'].to_numpy()) # [s]
            B_lambda.append(np.float64(df['Instrument Wavelength'].to_numpy()) * pow(10,-9)) # [m]

        return B_t, B_lambda

    def TC300(self, path):
        '''
            Read temperature measurements from TC300 data logs
        '''
        T_t, T_T1, T_T2 = [], [], []

        for file in sorted(path, key=self.sort_key):
            df = pd.read_csv(file, sep=',', header=None, skiprows=1, 
                             names=['Time', 'TargetTemp1', 'ActualTemp1', 'TargetCurrent1', 'ActualCurrent1', 'Voltage1', 'TargetTemp2', 'ActualTemp2', 'TargetCurrent2', 'ActualCurrent2', 'Voltage2'])
            
            df['Time'] = pd.to_datetime(df.iloc[:, 0], format='%H:%M:%S')
            df['Time'] = df['Time'].dt.hour * 3600 + df['Time'].dt.minute * 60 + df['Time'].dt.second

            T_t.append(df['Time'].to_numpy()) # [s]
            T_T1.append(df['ActualTemp1'].to_numpy()) # [°C]
            T_T2.append(df['ActualTemp2'].to_numpy()) # [°C]
        
        return T_t, T_T1, T_T2
    
    def Lock_ins(self, path):
        '''
            Read lock-ins data files
        '''
        para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod = [], [], [], [], [], [], [], []

        for file in sorted(path, key=self.sort_key):
            settings = []  # Store extracted values from the current file
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            value = parts[1]
                            # Skip 'Input' and 'gain'
                            if value.lower() not in ['input', 'gain']:
                                settings.append(float(value))  # Convert to float and append to list
                para.append(settings)
            
            data = pd.read_csv(file, sep=',', header=None, skiprows=8, 
                               names=['Timestamp', 'X_dc', 'Y_dc', 'X_2f', 'Y_2f', 'X_mod', 'Y_mod'])
            
            Lock_in_t.append(data['Timestamp'].to_numpy()) # [s]
            self.correct_negative_time_stamp(data['Timestamp'].to_numpy())
            Xdc.append(data['X_dc'].to_numpy()) # [V]
            Ydc.append(data['Y_dc'].to_numpy()) # [V]
            X2f.append(data['X_2f'].to_numpy()) # [V]
            Y2f.append(data['Y_2f'].to_numpy()) # [V]
            Xmod.append(data['X_mod'].to_numpy()) # [V]
            Ymod.append(data['Y_mod'].to_numpy()) # [V]
        
        para = np.array(para, dtype=object)
        return para, Lock_in_t, Xdc, Ydc, X2f, Y2f, Xmod, Ymod
    
    def comb_theta_vs_lambda(self, path):
        '''
            Read analyzed combined theta vs lambda data
        '''
        Date, N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste = [], [], [], [], [], [], [], [], []
        file_name = "Polarization_rotations_vs_wavelength_2.csv"
        file_path = os.path.join(path, file_name) 

        data = pd.read_csv(file_path, sep=',', header=None, skiprows=1, 
                           names=['Date', 'Number_of_Data_Points', 'Fitted_scan_amplitude', 'lambda_ave', 'lambda_std', 'lambda_ste', 'theta_ave', 'theta_std', 'theta_ste'])
        
        Date.extend(data['Date'])
        N.extend(data['Number_of_Data_Points'])
        A_fit.extend(data['Fitted_scan_amplitude'])
        Lambda_ave.extend(data['lambda_ave'])
        Lambda_std.extend(data['lambda_std'])
        Lambda_ste.extend(data['lambda_ste'])
        Theta_ave.extend(data['theta_ave'])
        Theta_std.extend(data['theta_std'])
        Theta_ste.extend(data['theta_ste'])
        
        return Date, N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste