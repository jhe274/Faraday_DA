import re
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
                             names=['Timestamp', 'Instrument Status',
                                    'Instrument Wavelength', 'Instrument Intensity'])
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
                             names=['Time', 'TargetTemp1', 'ActualTemp1',
                                    'TargetCurrent1', 'ActualCurrent1', 'Voltage1',
                                    'TargetTemp2', 'ActualTemp2', 'TargetCurrent2',
                                    'ActualCurrent2', 'Voltage2'])
            df['Time'] = pd.to_datetime(df.iloc[:, 0], format='%H:%M:%S')
            df['Time'] = df['Time'].dt.hour * 3600 + df['Time'].dt.minute * 60 + df['Time'].dt.second
            T_t.append(df['Time'].to_numpy()) # [s]
            T_T1.append(df['ActualTemp1'].to_numpy()) # [°C]
            T_T2.append(df['ActualTemp2'].to_numpy()) # [°C]
        
        return T_t, T_T1, T_T2
    
    def lockins(self, path):
        '''
        Read lock-ins data files
        '''
        para, lockins_t, Xmod, Ymod, X2f, Y2f, Xdc, Ydc = [], [], [], [], [], [], [], []

        for file in sorted(path, key=self.sort_key):
            settings = []                                                                   # Store extracted values from the current file
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            value = parts[1]
                            # Skip 'Input' and 'gain'
                            if value.lower() not in ['input', 'gain']:
                                settings.append(float(value))                               # Convert to float and append to list
                para.append(settings)
            df = pd.read_csv(file, sep=',', header=None, skiprows=9, 
                               names=['Timestamp', 'X_mod', 'Y_mod',
                                      'X_2f', 'Y_2f', 'X_dc', 'Y_dc'])
            df['Timestamp'] = pd.to_datetime(df['Timestamp'])
            start_time = df['Timestamp'].iloc[0]
            df['Timestamp'] = (df['Timestamp'] - start_time).dt.total_seconds()
            lockins_t.append(df['Timestamp'].to_numpy())                                    # [s]
            Xmod.append(df['X_mod'].to_numpy())                                             # [V]
            Ymod.append(df['Y_mod'].to_numpy())                                             # [V]
            X2f.append(df['X_2f'].to_numpy())                                               # [V]
            Y2f.append(df['Y_2f'].to_numpy())                                               # [V]
            Xdc.append(df['X_dc'].to_numpy())                                               # [V]
            Ydc.append(df['Y_dc'].to_numpy())                                               # [V]
        para = np.array(para, dtype=object)

        return para, lockins_t, Xmod, Ymod, X2f, Y2f, Xdc, Ydc