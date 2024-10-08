import re
import datetime as dt
import numpy as np
import pandas as pd

class Read:

    def sort_key(self, path):
        # Check if the file has an underscore and a numeric part at the end
        match = re.search(r'_(\d+)(?=\.[^.]*$|$)', path)
        if match:
            # Return the numeric part for files with underscore and number
            return int(match.group(1))
        
        # If no underscore and number, return a very low value to sort these files first
        return -1
        
    def DLCpro_WideScan(self, path):
        '''
        Read Toptica DLC pro wide scan output
        '''
        x, y, Y, DLCpro_t = [], [], [], []
        for file in sorted(path, key=self.sort_key):
            df = pd.read_csv(file, sep=',', header=None, skiprows=1,
                             names=['Piezo Voltage (V)', 'Fine In 1 (V)',
                                    'Monitor Photodiode (non-linear)', 'time (ms)'])
            x.append(df['Piezo Voltage (V)'].to_numpy())                                                                           # [V]
            y.append(df['Fine In 1 (V)'].to_numpy())                                                                               # [V]
            Y.append(df['Monitor Photodiode (non-linear)'].to_numpy())
            DLCpro_t.append(df['time (ms)'].to_numpy() * 1e-3)                                                                     # [s]
        
        return x, y, Y, DLCpro_t

    def Bristol(self, path):
        '''
        Read Wavelength measurements from Brilstol 871 wavelength meter
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
            B_t.append(df['Timestamp'].to_numpy())                                                                                 # [s]
            B_lambda.append(np.float64(df['Instrument Wavelength'].to_numpy()) * 1e-9)                                             # [m]
            
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
            T_t.append(df['Time'].to_numpy())                                                                                       # [s]
            T_T1.append(df['ActualTemp1'].to_numpy())                                                                               # [°C]
            T_T2.append(df['ActualTemp2'].to_numpy())                                                                               # [°C]
        
        return T_t, T_T1, T_T2

    def lockins(self, path):
        '''
        Read lock-ins data files
        '''
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = [], [], [], [], [], [], [], []
        for file in sorted(path, key=self.sort_key):
            settings = []                                                                                                           # Store extracted values from the current file
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            value = parts[1]
                            # Skip 'Input' and 'gain'
                            if value.lower() not in ['input', 'gain']:
                                settings.append(float(value))                                                                       # Convert to float and append to list
                para.append(settings)
            df = pd.read_csv(file, sep=',', header=None, skiprows=9, 
                               names=['Timestamp', 'X_1f', 'Y_1f',
                                      'X_2f', 'Y_2f', 'X_dc', 'Y_dc'])
            df['Timestamp'] = pd.to_datetime(df['Timestamp'])
            start_time = df['Timestamp'].iloc[0]
            df['Timestamp'] = (df['Timestamp'] - start_time).dt.total_seconds()
            lockins_t.append(df['Timestamp'].to_numpy())                                                                            # [s]
            X1f.append(df['X_1f'].to_numpy())                                                                                       # [V]
            Y1f.append(df['Y_1f'].to_numpy())                                                                                       # [V]
            X2f.append(df['X_2f'].to_numpy())                                                                                       # [V]
            Y2f.append(df['Y_2f'].to_numpy())                                                                                       # [V]
            Xdc.append(df['X_dc'].to_numpy())                                                                                       # [V]
            Ydc.append(df['Y_dc'].to_numpy())                                                                                       # [V]
        para = np.array(para, dtype=object)

        return para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc
    
    def ellip_theta(self, path):
        '''
        Read wavelength, ellipticities and Faraday rotations
        '''
        def read_header(file):
            df = pd.read_csv(file, header=None, usecols=[1], nrows=4)
            return df.iloc[0, 0], df.iloc[1, 0], df.iloc[2, 0], df.iloc[3, 0]

        def read_data(file):
            df = pd.read_csv(file, sep=',', header=None, skiprows=5,
                            names=['Wavelength (m)', 'Ellipticity (radian)', 'Faraday rotation (radian)'])
            wl = df['Wavelength (m)'].to_numpy(dtype=np.float64)
            ellip = df['Ellipticity (radian)'].to_numpy(dtype=np.float64)
            theta = df['Faraday rotation (radian)'].to_numpy(dtype=np.float64)
            return wl, ellip, theta

        date, temp, Bz, power, wl, ellip, theta = [], [], [], [], [], [], []

        for file in sorted(path, key=self.sort_key):
            d, t, b, p = read_header(file)
            date.append(d)                                                                                                          # [XX-XX-XXXX]
            temp.append(t)                                                                                                          # [°C]
            Bz.append(b)                                                                                                            # [G]
            power.append(p)                                                                                                         # [microW]

            w, e, th = read_data(file)
            wl.append(w)                                                                                                            # [m]
            ellip.append(e)                                                                                                         # [rad]
            theta.append(th)                                                                                                        # [rad]

        return date, temp, Bz, power, wl, ellip, theta