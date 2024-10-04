import os, glob
import datetime as dt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze
from scipy.signal import find_peaks

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()
    
    def find_closest_indices(self, arrays, power, n):
        arrays = [np.asarray(subarray, dtype=np.float64) for subarray in arrays]
        diffs = [np.abs(subarray - power) for subarray in arrays]
        
        # Create a structured array to keep track of subarray index, element index, value, and difference
        all_closest = np.array([(sub_idx, idx, subarray[idx], diff[idx])
                                for sub_idx, (subarray, diff) in enumerate(zip(arrays, diffs))
                                for idx in range(len(subarray))],
                            dtype=[('sub_idx', int), ('idx', int), ('value', float), ('diff', float)])
        
        # Sort by the difference and select the top `n` closest values
        sorted_closest = np.sort(all_closest, order='diff')[:n]
        
        # Extract the result as a list of tuples
        result = [(item['sub_idx'], item['idx'], item['value']) for item in sorted_closest]
        return result
    
    def corr_plot(self, dates, phytype, datatype, power, n):
        
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        results = [self.reader.ellip_theta(glob.glob(os.path.join(processed_path, date, '*.csv'))) for date in dates]
        T, B, P, x1, y1, y2 = zip(*[(res[1], res[2], res[3], res[4], res[5], res[6]) for res in results])

        # Convert to float arrays
        T, B, P = map(self.analyzer.convert_to_float_array, [T, B, P])
        
        # Find closest indices
        result = self.find_closest_indices(P, power, n)
    
        for sub_idx, idx, value in result:
            nu = self.consts.c / x1[sub_idx][idx]                                                                                       # [Hz]
            detun = (nu - self.consts.Nu39_D2) * 1e-9                                                                                   # [GHz]
            smoothed_y1 = self.analyzer.smooth(y1[sub_idx][idx]*1e6, 90)
            smoothed_y2 = self.analyzer.smooth(y2[sub_idx][idx]*1e6, 90)
            y1_grad = np.gradient(y1[sub_idx][idx]*1e6, detun)
            smoothed_y1_grad = np.gradient(smoothed_y1, detun)
            y2_grad = np.gradient(y2[sub_idx][idx]*1e6, detun)
            smoothed_y2_grad = np.gradient(smoothed_y2, detun)
            label = f'T={T[sub_idx][idx]:.2f}°C, $B_z$={B[sub_idx][idx]:.2f} G, P={P[sub_idx][idx]:.2f} μW'

            if phytype == 'CD':
                if datatype == 'raw':
                    ax.plot(detun, y1[sub_idx][idx] * 1e6, '.', label=label, markersize=4)
                    # ax.plot(detun, smoothed_y1, '.', label=label, markersize=4)
                elif datatype == 'grad':
                    ax.plot(detun, y1_grad, '.', label=label, markersize=4)
                    # ax.plot(detun, smoothed_y1_grad, '.', label=label, markersize=4)
            elif phytype == 'CB':
                if datatype == 'raw':
                    # peaks, _ = find_peaks(y2[sub_idx][idx] * 1e6, prominence=900, height=(2000, 3000))
                    # valleys, _ = find_peaks(-y2[sub_idx][idx] * 1e6, prominence=500, height=(2000, 2500))
                    # print(peaks, valleys)
                    # for index, indices in enumerate([peaks, valleys]):
                    #     for k in indices:
                    ax.plot(detun, y2[sub_idx][idx] * 1e6, '.', label=label, markersize=4)
                            # ax.axvline(x=detun[k], linestyle='--', label=label, markersize=4)
                            # print(detun[k])
        
                    if sub_idx == 0:
                        initial_guess = [0.8, 21.85, -4.05, 0, 45*1e-6]
                        bounds = ([.1, 21.7, -4.15, -1, -100*1e-6], [2, 23, -3.95, 1, 100*1e-6])
                        params, covariance = curve_fit(self.theory.resonant_FR, nu, y2[sub_idx][idx], p0=initial_guess, bounds=bounds)
                        Kn_fit, T_fit, B_fit, P_fit, const = [np.round(val,4) for val in params[:4]] + [np.round(params[4] * 1e6,4)]
                        print(Kn_fit, T_fit, B_fit, P_fit, const)
                        label_fit = f'[K]={np.round(Kn_fit,2)}$\\times10^{{14}}$ m$^{{-3}}$, $B_z$={np.round(B_fit,2)} G, $P$={np.round(P_fit*1e2,2)}%'
                        ax.plot(detun, self.theory.resonant_FR(nu, Kn_fit, T_fit, B_fit, P_fit, const*1e-6)*1e6, '--', 
                                label=label_fit, linewidth=3)
                        
                    elif sub_idx == 14:
                        initial_guess = [1.4, 20.3, -6.11, 0, 30*1e-6]
                        bounds = ([.05, 20.3, -6.3, -1, -100*1e-6], [5, 23, -5.9, 1, 100*1e-6])
                        params, covariance = curve_fit(self.theory.resonant_FR, nu, y2[sub_idx][idx], p0=initial_guess, bounds=bounds)
                        Kn_fit, T_fit, B_fit, P_fit, const = [np.round(val,4) for val in params[:4]] + [np.round(params[4] * 1e6,4)]
                        print(Kn_fit, T_fit, B_fit, P_fit, const)
                        label_fit = f'[K]={np.round(Kn_fit,2)}$\\times10^{{14}}$ m$^{{-3}}$, $B_z$={np.round(B_fit,2)} G, $P$={np.round(P_fit*1e2,2)}%'
                        ax.plot(detun, self.theory.resonant_FR(nu, Kn_fit, T_fit, B_fit, P_fit, const*1e-6)*1e6, '--', 
                                label=label_fit, linewidth=3)
                        
                    elif sub_idx == 6:
                        initial_guess = [1.1, 20.7, -5.12, 0, 34*1e-6]
                        bounds = ([.5, 20.7, -5.3, -1, -100*1e-6], [1.5, 23, -4.9, 1, 100*1e-6])
                        params, covariance = curve_fit(self.theory.resonant_FR, nu, y2[sub_idx][idx], p0=initial_guess, bounds=bounds)
                        Kn_fit, T_fit, B_fit, P_fit, const = [np.round(val,4) for val in params[:4]] + [np.round(params[4] * 1e6,4)]
                        print(Kn_fit, T_fit, B_fit, P_fit, const)
                        label_fit = f'[K]={np.round(Kn_fit,2)}$\\times10^{{14}}$ m$^{{-3}}$, $B_z$={np.round(B_fit,2)} G, $P$={np.round(P_fit*1e2,2)}%'
                        ax.plot(detun, self.theory.resonant_FR(nu, Kn_fit, T_fit, B_fit, P_fit, const*1e-6)*1e6, '--', 
                                label=label_fit, linewidth=3)
                    # ax.plot(detun, self.theory.resonant_FR(nu, 0.83, 21.85, -4.05, 0.013, 43.2*1e-6)*1e6, '--', color='red', label='Manual fit')
                    # ax.plot(detun, self.theory.resonant_FR(nu, 1.474, 19.497, -5.12, .002, 60*1e-6)*1e6, '--', color='red', label='Manual fit')
                elif datatype == 'grad':
                    ax.plot(detun, y2_grad, '.', label=label, markersize=4)
                    # ax.plot(detun, smoothed_y2_grad, '.', label=label, markersize=4)

                    peaks, _ = find_peaks(smoothed_y2_grad, prominence=600, height=(4500, 5100))
                    valleys, _ = find_peaks(-smoothed_y2_grad, prominence=1000, height=(3500, 4000))
                    print(peaks, valleys)

                    for idx, indices in enumerate([peaks, valleys]):
                        for k in indices:
                            ax.axvline(x=detun[k], linestyle='--')
                            print(detun[k])

        self.plot_settings(power, phytype, datatype)
    
    def plot_settings(self, power, phytype, datatype):
        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.yticks(fontsize=25)
        plt.xlim(-5,6)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        plt.legend(loc='best', markerscale=10, fontsize=25)

        plots_path = os.path.join(processed_path, 'Correlation plots', 'Sort by power', f'{power:.0f} Microwatts')
        os.makedirs(plots_path, exist_ok=True)

        if phytype == 'CD':
            if datatype == 'raw':
                plt.ylabel(r'Ellipticity (microrad.)', fontsize=25)
                plt.title(f'Ellipticity vs Frequency Detuning with $P\simeq${power:.0f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_Ellipticity_vs_Detuning_{power:.0f}microW.png'))
            elif datatype == 'grad':
                plt.ylabel(r'$\nabla_{\nu}\epsilon$ (microrad/GHz)', fontsize=25)
                plt.title(rf'$\nabla_\nu\epsilon$ vs Frequency Detuning with $P\simeq${power:.0f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_gradEllipticity_vs_Detuning_{power:.0f}microW.png'))
        elif phytype == 'CB':
            if datatype == 'raw':
                plt.ylabel(r'Faraday Rotation (microrad.)', fontsize=25)
                plt.title(f'Faraday Rotation vs Frequency Detuning with $P\simeq${power:.2f} μW', fontsize=25)
                # plt.savefig(os.path.join(plots_path, f'[X]_FR_vs_Detuning_{power:.0f}microW(fitted).png'))
            elif datatype == 'grad':
                plt.ylabel(r'$\nabla_{\nu}\theta$ (microrad/GHz)', fontsize=25)
                plt.title(rf'$\nabla_\nu\theta$ vs Frequency Detuning with $P\simeq${power:.0f} μW', fontsize=25)
                plt.savefig(os.path.join(plots_path, f'[X]_gradFR_vs_Detuning_{power:.0f}microW(smoothed).png'))
        plt.show()

if __name__ == "__main__":
    dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    # dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed data')

    plotter = Plot()
    dates = ['06-18-2024', '06-24-2024', '06-25-2024', '06-26-2024', '06-27-2024', '07-01-2024',
                '05-07-2024', '05-09-2024', '05-15-2024', '05-19-2024',
                '05-23-2024', '05-29-2024', '05-31-2024', '06-05-2024', '06-07-2024']
    plotter.corr_plot(dates, 'CB', 'raw', 0.5, 3)