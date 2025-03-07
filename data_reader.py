import re
import numpy as np
import pandas as pd

class DataReader:
    """
    A class to read and process various experimental data files related to 
    laser scans, wavelength measurements, temperature logs, lock-in amplifier 
    signals, and processed data analysis.
    """

    @staticmethod
    def sort_key(path):
        """
        Extracts the numeric part from filenames to ensure proper sorting.
        Files with numeric suffixes are sorted numerically; others appear first.
        """
        match = re.search(r'_(\d+)(?=\.[^.]*$|$)', path)
        return int(match.group(1)) if match else -1
    
    @staticmethod
    def read_csv_data(path, **kwargs):
        """
        Reads multiple CSV files, sorts them using `sort_key`, and loads them into DataFrames.
        Returns a list of DataFrames.
        """
        return [pd.read_csv(file, **kwargs) for file in sorted(path, key=DataReader.sort_key)]
        
    def read_dlcpro_widescan(self, path):
        """
        Reads and processes Toptica DLC Pro wide scan output files.
        Returns: 
            - x: Piezo voltage array [V]
            - y1: Fine In 1 array [V], depending on the settings
            - y2: Monitor photodiode array (non-linear), depending on the settings
            - timestamp: Time array [s]
        """
        dfs = self.read_csv_data(path, sep=',', header=None, skiprows=1,
                                 names=['Piezo Voltage (V)', 'Fine In 1 (V)', 
                                        'Monitor Photodiode (non-linear)', 'time (ms)'])
        x, y1, y2, timestamp = [], [], [], []
        for df in dfs:
            x.append(df['Piezo Voltage (V)'].to_numpy())
            y1.append(df['Fine In 1 (V)'].to_numpy())
            y2.append(df['Monitor Photodiode (non-linear)'].to_numpy())
            timestamp.append(df['time (ms)'].to_numpy() * 1e-3)

        return x, y1, y2, timestamp
    
    def read_bristol(self, path):
        """
        Reads wavelength measurements from Bristol 871 wavelength meter.
        Returns:
            - timestamp: Time array [s]
            - wavelength: Wavelength array [m]
        """
        dfs = self.read_csv_data(path, sep=',', header=None, skiprows=1,
                                 names=['Timestamp', 'Instrument Status', 
                                        'Instrument Wavelength', 'Instrument Intensity'])
        timestamp, wavelength = [], []
        for df in dfs:
            t0 = pd.to_datetime(df['Timestamp']).iloc[0]  # Reference timestamp
            timestamp.append(pd.to_datetime(df['Timestamp']).sub(t0).dt.total_seconds().to_numpy())
            wavelength.append(df['Instrument Wavelength'].to_numpy(dtype=np.float64) * 1e-9)

        return timestamp, wavelength
    
    def read_gaussmeter(self, path):
        """
        Reads magnetic field measurements from Lakeshore 475 DSP Gaussmeter.
        Returns:
            - timestamp: Time array [s]
            - B0s: Magnetic field components [G]
            - temps: Temperature array [°C]
        """
        dfs = self.read_csv_data(path, sep=',', header=None, skiprows=1,
                                 names=['Timestamp','MagneticFluxDensity(G)','Temperature(C)'])
        timestamps, B0s, temps = [], [], []
        for df in dfs:
            t0 = pd.to_datetime(df['Timestamp']).iloc[0]
            timestamps.append(pd.to_datetime(df['Timestamp']).sub(t0).dt.total_seconds().to_numpy())
            B0s.append(df['MagneticFluxDensity(G)'].to_numpy())
            temps.append(df['Temperature(C)'].to_numpy())
        
        return timestamps, B0s, temps

    def read_tc300(self, path):
        """
        Reads temperature measurements from TC300 data logs.
        Returns:
            - timestamp: Time array [s]
            - temp_T1: Actual temperature sensor 1 [°C]
            - temp_T2: Actual temperature sensor 2 [°C]
        """
        dfs = self.read_csv_data(path, sep=',', header=None, skiprows=1,
                                 names=['Time', 'TargetTemp1', 'ActualTemp1', 'TargetCurrent1', 'ActualCurrent1', 'Voltage1',
                                        'TargetTemp2', 'ActualTemp2', 'TargetCurrent2', 'ActualCurrent2', 'Voltage2'])
        timestamp, temp_T1, temp_T2 = [], [], []
        for df in dfs:
            time_dt = pd.to_datetime(df['Time'], format='%H:%M:%S')
            time_seconds = time_dt.dt.hour * 3600 + time_dt.dt.minute * 60 + time_dt.dt.second
            timestamp.append(time_seconds.to_numpy())
            temp_T1.append(df['ActualTemp1'].to_numpy())
            temp_T2.append(df['ActualTemp2'].to_numpy())

        return timestamp, temp_T1, temp_T2

    def read_lockins(self, path):
        """
        Reads lock-in amplifier data files.
        Returns:
            - para: Extracted metadata parameters
            - timestamp: Time array [s]
            - X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f: Lock-in measurement arrays [V]
        """
        para, timestamp, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f = [], [], [], [], [], [], [], [], [], []
        
        for file in sorted(path, key=self.sort_key):
            settings = []
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        parts = line.strip('#').strip()
                        
                        # Old format: "#SENS_1f[V] 0.01"
                        if re.match(r'.*\s[\d\.Ee+-]+$', parts):
                            key, value = parts.rsplit(' ', 1)
                        # New format: "#1f Sensitivity [V]: 0.002"
                        elif ':' in parts:
                            key, value = parts.split(':', 1)
                        else:
                            continue  # Skip lines without numerical metadata
                        
                        try:
                            settings.append(float(value.strip()))
                        except ValueError:
                            pass  # Skip non-numeric metadata entries
                    
            para.append(settings)

            # Determine correct skiprows value dynamically
            with open(file, 'r') as f:
                lines = f.readlines()
                skiprows = sum(1 for line in lines if line.startswith('#')) + 1

            df = pd.read_csv(file, sep=',', header=None, skiprows=skiprows, 
                            names=['Timestamp', 'X_1f', 'Y_1f', 'X_2f', 'Y_2f', 'X_dc', 'Y_dc', 'X_m2f', 'Y_m2f'])
            
            df['Timestamp'] = pd.to_datetime(df['Timestamp'])
            df['Timestamp'] = (df['Timestamp'] - df['Timestamp'].iloc[0]).dt.total_seconds()
            
            # Compute step size for downsampling
            step_size = int(5)  # take every 5th point
            
            # Apply downsampling
            timestamp.append(df['Timestamp'].iloc[::step_size].to_numpy())
            X1f.append(df['X_1f'].iloc[::step_size].to_numpy())
            Y1f.append(df['Y_1f'].iloc[::step_size].to_numpy())
            X2f.append(df['X_2f'].iloc[::step_size].to_numpy())
            Y2f.append(df['Y_2f'].iloc[::step_size].to_numpy())
            Xdc.append(df['X_dc'].iloc[::step_size].to_numpy())
            Ydc.append(df['Y_dc'].iloc[::step_size].to_numpy())
            Xm2f.append(df['X_m2f'].iloc[::step_size].to_numpy())
            Ym2f.append(df['Y_m2f'].iloc[::step_size].to_numpy())

        return para, timestamp, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f
    
    def read_processed_da(self, path):
        """
        Reads processed data analysis files containing wavelength, ellipticity, and Faraday rotation.
        Returns:
            - date: Date of the experiment
            - temp: Temperature [°C]
            - Bz: Longitudinal magnetic field strength [G]
            - power: Laser power [µW]
            - wl: Wavelength array [m]
            - ellipticity: Ellipticity array [rad]
            - angle: Faraday rotation angle array [rad]
        """
        def read_header(file):
            df = pd.read_csv(file, header=None, usecols=[1], nrows=4)
            return df.iloc[:, 0].tolist()

        def read_data(file):
            df = pd.read_csv(file, sep=',', header=None, skiprows=5,
                             names=['Wavelength (m)', 'Ellipticity (radian)', 'Faraday rotation (radian)'])
            return df['Wavelength (m)'].to_numpy(), df['Ellipticity (radian)'].to_numpy(), df['Faraday rotation (radian)'].to_numpy()

        date, temp, Bz, power, wl, ellip, angle = [], [], [], [], [], [], []

        for file in sorted(path, key=self.sort_key):
            d, t, b, p = read_header(file)
            date.append(d)
            temp.append(t)
            Bz.append(b)
            power.append(p)

            w, e, a = read_data(file)
            wl.append(w)
            ellip.append(e)
            angle.append(a)

        return date, temp, Bz, power, wl, ellip, angle