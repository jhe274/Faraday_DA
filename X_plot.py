import os, glob
import datetime as dt
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()

    def read_data(self, lambda_path, lockin_path, run, n):
        """
        Read measured wavelength and voltages and calculate the ellipticities and Faraday rotations
        """
        B_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockin_path)
        epsilon, theta = self.analyzer.FR_double_Kvapor(lockin_path, X1f, X2f, Xdc)
        x, x0, Eps, The = [], [], [], []

        runs = range(run-1, run+1)

        for i in runs:
            B_t[i], Lambda[i] = self.analyzer.filter_data(B_t[i], Lambda[i])

            B_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], epsilon[i])
            B_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], theta[i])

            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[i], lockins_t[i], para[i][2], n)
            Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])

            x0.append(Lambd)                                                                                                                                # [m]
            x.append(self.consts.c / x0[i - run + 1] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                   # [GHz]
            Eps.append(ep * 1e3)                                                                                                                            # [millirad]
            The.append(th * 1e3)                                                                                                                            # [microrad]
            
        return x0, x, Eps, The
    
    def physics_extraction(self, lambda_path, lockin_path, run, n):
        """
        Calculate ellipticity/Faraday rotation by subtraction appropriate background
        """
        x0, x, Eps, The = self.read_data(lambda_path, lockin_path, run, n)
        CD_empty, CD_vapor, CD_K, CB_empty, CB_vapor, CB_K = [], [], [], [], [], []

        runs = range(run-1, run+1)

        if len(runs) > 2:
            for idx, x_val in enumerate(x0[0]):
                # Considr air ellipticity/FR as background measurement
                idx_empty = np.argmin(np.abs(x0[1] - x_val))
                idx_vapor = np.argmin(np.abs(x0[2] - x_val))
                CD_empty.append(Eps[1][idx_empty] - Eps[0][idx])
                CD_vapor.append(Eps[2][idx_vapor] - Eps[0][idx])
                CB_empty.append(The[1][idx_empty] - The[0][idx])
                CB_vapor.append(The[2][idx_vapor] - The[0][idx])
            for idx, x_val in enumerate(x0[1]):
                # Consider empty cell ellipticity/FR as background measurement
                idx_K = np.argmin(np.abs(x0[2] - x_val))
                CD_K.append(Eps[2][idx_K] - Eps[1][idx])
                CB_K.append(The[2][idx_K] - The[1][idx])
        else:
            for idx, x_val in enumerate(x0[0]):
                # Consider air ellipticity/FR as background measurement
                idx_vapor = np.argmin(np.abs(x0[1] - x_val))
                CD_vapor.append(Eps[1][idx_vapor] - Eps[0][idx])
                CB_vapor.append(The[1][idx_vapor] - The[0][idx])

        return x, CD_empty, CB_empty, CD_vapor, CB_vapor, CD_K, CB_K
    
    def raw_plot(self, lambda_path, lockin_path, run, n, B, power, dtype):
        """
        Plot ellipticity/FR for each individual measurements
        """
        x0, x, Eps, The = self.read_data(lambda_path, lockin_path, run, n)
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        runs = range(run-1, run+1)

        if len(runs) > 2:
            for i in runs:
                if dtype == 'CD':
                    ax.plot(x[i - run + 1], Eps[i - run + 1][1:], label=(r'$\epsilon_\text{air}$' if i-run+1 == 0 
                                                                else (r'$\epsilon_\text{empty cell}$' if i-run == 0
                                                                    else r'$\epsilon_\text{vapor cell}$')
                                                                )
                                )
                elif dtype == 'CB':
                    ax.plot(x[i - run + 1], The[i - run + 1][1:], label=(r'$\theta_\text{air}$' if i-run+1 == 0 
                                                            else (r'$\theta_\text{empty cell}$' if i-run == 0
                                                                    else r'$\theta_\text{vapor cell}$')
                                                            )
                            )
        else:
            for i in runs:
                if dtype == 'CD':
                    ax.plot(x[i - run + 1], Eps[i - run + 1][1:], label=(r'$\epsilon_\text{air}$' if i-run+1 == 0 
                                                        else r'$\epsilon_\text{vapor cell}$')
                        )
                elif dtype == 'CB':
                    ax.plot(x[i - run + 1], The[i - run + 1][1:], label=(r'$\theta_\text{air}$' if i-run+1 == 0 
                                                            else r'$\theta_\text{vapor cell}$')
                            )
        
        self.plot_settings(n, B, power, date, dtype)

    def extracted_plot(self, lambda_path, lockin_path, run, n, B, power, dtype, material):
        """
        Plot background subtracted ellipticity/FR
        """
        x, CD_empty, CB_empty, CD_vapor, CB_vapor, CD_K, CB_K = self.physics_extraction(lambda_path, lockin_path, run, n)
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        if dtype == 'CD':
            if material == 'empty':
                ax.plot(x[0], CD_empty, '.', label=r'$\epsilon_\text{empty cell}-\epsilon_\text{air}$', markersize=2)
            elif material == 'vapor':
                ax.plot(x[0], CD_vapor, '.', label=r'$\epsilon_\text{vapor cell}-\epsilon_\text{air}$', markersize=2)
            elif material == 'K':
                ax.plot(x[1], CD_K, '.', label=r'$\epsilon_\text{vapor cell}-\epsilon_\text{empty cell}$', markersize=2)
        elif dtype == 'CB':
            if material == 'empty':
                ax.plot(x[0], CB_empty, '.', label=r'$\theta_\text{empty cell}-\theta_\text{air}$', markersize=2)
            elif material == 'vapor':
                ax.plot(x[0], CB_vapor, '.', label=r'$\theta_\text{vapor cell}-\theta_\text{air}$', markersize=2)
            elif material == 'K':
                ax.plot(x[1], CB_K, '.', label=r'$\theta_\text{vapor cell}-\theta_\text{empty cell}$', markersize=2)

        # self.peaks_valleys_plot(x[0], CB_vapor)
        self.plot_settings(run, B, power, date, dtype)

    def peaks_valleys_plot(self, x, y):
        """
        Find peaks and valleys in the ellipticity/FR
        """
        peaks, _ = find_peaks(np.array(y), prominence=20, height=(200, 600))
        valleys, _ = find_peaks(-1*np.array(y), prominence=10, height=(-600, -200))

        for idx, indices in enumerate([peaks, valleys]):
            for k in indices:
                plt.vlines(x=x[k], ymin=0, ymax=y[k], linestyle=':', color='purple')

    def curve_fitting(self):
        l = (7.5-0.159*2)*1e-2                                                                                                          # [m]
        nu_D1 = 389286.058716 * 1e9
        nu_D2 = 391016.17003 * 1e9

        def FR(nu, Kn, T, B, P, const):
            delta_nu_D2 = nu - nu_D2
            delta_nu_D1 = nu - nu_D1
            delta_doppler_D2 = self.theory.doppler_broad(nu_D2, T)
            delta_doppler_D1 = self.theory.doppler_broad(nu_D1, T)
            term1 = (
                (7*(delta_nu_D2**2 - delta_doppler_D2**2/4) / (delta_nu_D2**2 + delta_doppler_D2**2/4)**2) + 
                (4*(delta_nu_D1**2 - delta_doppler_D1**2/4) / (delta_nu_D1**2 + delta_doppler_D1**2/4)**2) - 
                (2*(delta_nu_D1*delta_nu_D2) / ((delta_nu_D2-delta_doppler_D2) * (delta_nu_D1-delta_doppler_D1))**2)) / (3 * self.consts.h)
            # term2 = np.sign(B) * ((nu / (nu_D1 * (nu - nu_D1))) - (nu / (nu_D2 * (nu - nu_D2)))) / (self.consts.k_B * T)
            
            diamagnetic_theta = self.consts.mu_B * B*1e-4 * (term1)
            paramagnetic_thetea = P * (
                    (delta_nu_D2 / ((delta_nu_D2 - delta_doppler_D2) ** 2)) -
                    (delta_nu_D1 / ((delta_nu_D1 - delta_doppler_D1) ** 2))
                )
            # paramagnetic_thetea = P * (
            #         (delta_nu_D2 / ((delta_nu_D2 - doppler_broad_D2)**2 + (doppler_broad_D2**2)/4)) -
            #         (delta_nu_D1 / ((delta_nu_D1 - doppler_broad_D1)**2 + (doppler_broad_D1**2)/4))
            #     )
            theta = self.consts.alpha * Kn*1e14 * l * (diamagnetic_theta + paramagnetic_thetea) * 1e6 + const                                                        # [microrad]
            return theta

        initial_guess = [1.47, 19.5, -5.103, -0.002, -60]
        bounds = ([.1, 15, -5.2, -1, -100], [5, 25, -5., 1, 100])
        # params, covariance = curve_fit(FR, x, y, p0=initial_guess, bounds=bounds)
        # print(params)
        # Kn, T, Bz, PK, const = np.round(params,3)
        # plt.plot(x*1e-9 - nu_D2 * 1e-9, FR(x, Kn, T, Bz, PK, const), '--', color='red', label='Curve fit')
        # plt.plot(x*1e-9 - nu_D2 * 1e-9, FR(x, 1.474, 19.497, -5.103, -.002, -60), '--', color='green', label='Manual fit')

    def plot_settings(self, run, B, power, date, dtype):
        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # plt.ylim(400,-650)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        plt.legend(loc='best', fontsize=25)
        if dtype == 'CD':
            plt.ylabel(r'Ellipticity (millirad.)', fontsize=25)
            plt.title(f'Ellipticity vs Frequency, $B_z$={-B} G, $P$={power} $\mu$W @{date}', fontsize=25)
            plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity_vs_Frequency_{date}_run{run}-{run+1}.png'))
        elif dtype == 'CB':
            plt.ylabel(r'Faraday Rotation (millirad.)', fontsize=25)
            plt.title(f'Faraday Rotation vs Frequency, $B_z$={-B} G, $P$={power} $\mu$W @{date}', fontsize=25)
            plt.savefig(os.path.join(Plots, f'{date}', f'FR_vs_Frequency_{date}_run{run}-{run+1}.png'))
        # plt.title(rf'$n={Kn}\times10^{{14}}\text{{m}}^3$, $T={T}^\circ$C, $B_z={Bz}$G, $P=.2\%$, $\theta_\text{{offset}}={const}\mu\text{{rad}}$', fontsize=25)
        plt.show()

if __name__ == "__main__":
    # dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
    K_vapor = os.path.join(dir_path, 'K vapor cell')
    Bristol = os.path.join(K_vapor, 'Bristol data')
    Lockins = os.path.join(K_vapor, 'Lockins data')
    Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '05-24-2024'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
    Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
    plotter.extracted_plot(Bristol_path, Lockins_path, 1, 5, 6.07, 99.8, 'CD', 'vapor')