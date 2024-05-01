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

dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
# dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
K_vapor = os.path.join(dir_path, 'K vapor cell')
Bristol = os.path.join(K_vapor, 'Bristol data')
Lockins = os.path.join(K_vapor, 'Lockins data')
Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()

    def Ellipticity_FR_X(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        B_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockin_path)
        epsilon, theta = self.analyzer.FR_double_Kvapor(lockin_path, X1f, X2f, Xdc)
        
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        x, x0, Eps, The = [], [], [], []

        """
        Plot measured epsilon/theta for vapor cell and background
        """
        for i in range(run-1, run+1):
            B_t[i], Lambda[i] = self.analyzer.filter_data(B_t[i], Lambda[i])

            B_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], epsilon[i])
            B_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], theta[i])

            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[i], lockins_t[i], para[i][2], n)
            Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])

            x0.append(Lambd)
            Eps.append(ep * 1e3)
            The.append(th * 1e6)
            # x = self.consts.c / x0[i - run + 1] * 1e-9 - self.consts.Nu39_D2 * 1e-9                                                                   # [GHz]

            # colors = 'y' if i == run - 1 else ('r' if i == run else 'b')
            # ax.plot(x, Eps[i - run + 1][1:], color=colors, label=(r'$\epsilon_\text{empty cell}$' if i-run+1 == 0 
            #                                         else (r'$\epsilon_\text{vapor cell}$' if i-run == 0
            #                                               else r'$\epsilon_\text{air}$')
            #                                         )
            #         )
            # ax.plot(x, The[i - run + 1][1:], color=colors, label=(r'$\theta_\text{empty cell}$' if i-run+1 == 0 
            #                                         else (r'$\theta_\text{vapor cell}$' if i-run == 0
            #                                               else r'$\theta_\text{air}$')
            #                                         )
            #         )

        """
        Subtracting measured background epsilon/theta values from sample epsilon/theta values
        """
        diff_eps, diff_the = [], []
        short_idx, long_idx = (0, 1) if len(x0[0]) <= len(x0[1]) else (1, 0)
        
        for idx, x_val in enumerate(x0[short_idx]):
            idx_long = np.argmin(np.abs(x0[long_idx] - x_val))
            diff_eps.append(Eps[short_idx][idx] - Eps[long_idx][idx_long])
            # diff_eps.append(Eps[long_idx][idx_long] - Eps[short_idx][idx])
            diff_the.append(The[short_idx][idx] - The[long_idx][idx_long])
            # diff_the.append(The[long_idx][idx_long] - The[short_idx][idx])

        # x0.append(self.consts.c / x0[short_idx] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                            # [GHz]
        # ax.plot(x0[2], diff_eps, label=r'$\epsilon_\text{K}$')
        # ax.plot(x0[2], diff_the, label=r'$\theta_{K}$')

        x = self.consts.c / x0[short_idx]
        y = np.array(diff_the)
        ax.plot(x*1e-9-self.consts.Nu39_D2*1e-9, y, '.', label=r'Measured $\theta_{K}$', color='blue', markersize=1)

        """
        Find the peaks and valleys in the epsilon_K/theta_K and mirror the data w.r.t. the peak axis
        """
        # peaks, _ = find_peaks(y, prominence=20, height=(-150, -100))
        # valleys, _ = find_peaks(-y, prominence=50, height=(100, 400))

        # for idx, indices in enumerate([peaks, valleys]):
        #     for k in indices:
        #         y0 = y[k]
        #         color = 'cyan' if idx == 0 else 'purple'
        #         ax.vlines(x=x[k]*1e-9-self.consts.Nu39_D2*1e-9, ymin=-400, ymax=y[k], linestyle=':', color=color)
        
        # mirror = -y[peaks[0]:] + 2 * y[peaks[0]]
        # new_y = y.copy()
        # new_y[peaks[0]:] = mirror
        # ax.plot(x[peaks[0]:]*1e-9-self.consts.Nu39_D2*1e-9, mirror, label=r'Manipulated $\theta_{K}$', color='green')

        """
        Curve fitting
        """
        def Kn_density(T):
            Kp = 10 ** (9.967 - 4646 / (273.15 + T))
            return Kp / (self.consts.k_B * (273.15 + T))
        # print(Kn_density(22)) # 4.13 * 1e14
        # print(Kn_density(24)) # 5.23 * 1e14
        # print(Kn_density(25)) # 5.88 * 1e14
        # print(Kn_density(26)) # 6.61 * 1e14

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
        params, covariance = curve_fit(FR, x, y, p0=initial_guess, bounds=bounds)
        print(params)
        Kn, T, Bz, PK, const = np.round(params,3)
        # plt.plot(x*1e-9 - nu_D2 * 1e-9, FR(x, Kn, T, Bz, PK, const), '--', color='red', label='Curve fit')
        plt.plot(x*1e-9 - nu_D2 * 1e-9, FR(x, 1.474, 19.497, -5.103, -.02, -60), '--', color='green', label='Manual fit')

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        # plt.ylabel(r'Ellipticity (millirad.)', fontsize=25)
        plt.ylabel(r'Faraday Rotation (microrad.)', fontsize=25)
        # plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        plt.ylim(150,-600)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        # plt.title(f'Ellipticity vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} $\mu$W @{date}', fontsize=25)
        # plt.title(f'Faraday Rotation vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.title(rf'$n={Kn}\times10^{{14}}\text{{m}}^3$, $T={T}^\circ$C, $B_z={Bz}$G, $P=.2\%$, $\theta_\text{{offset}}={const}\mu\text{{rad}}$', fontsize=25)
        # plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity(X)_vs_Frequency_{date}_run{run}-{run+1}.png'))
        plt.savefig(os.path.join(Plots, f'{date}', f'FR_vs_Frequency_{date}_run{run}-{run+1}(fitted).png'))
        plt.show()
    
    def Ellipticity_FR_R(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        B_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc = self.analyzer.R_lockins(lockin_path)
        epsilon, theta = self.analyzer.FR_double_Kvapor(lockin_path, R1f, R2f, Rdc)
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        x, x0, R1, R2, R3, Eps, The = [], [], [], [], [], [], []

        for i in range(run-1, run+1):
            B_t[i], Lambda[i] = self.analyzer.filter_data(B_t[i], Lambda[i])

            # trim R values
            B_t[i], Lambda[i], lockins_t[i], R1f[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], R1f[i])
            B_t[i], Lambda[i], lockins_t[i], R2f[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], R2f[i])
            B_t[i], Lambda[i], lockins_t[i], Rdc[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], Rdc[i])
            B_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], epsilon[i])
            B_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], theta[i])

            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[i], lockins_t[i], para[i][2], n)

            # average R values
            Lambd, r1 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R1f[i][l_idx])
            Lambd, r2 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R2f[i][l_idx])
            Lambd, r3 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], Rdc[i][l_idx])
            Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])

            x0.append(Lambd)
            R1.append(r1)
            R2.append(r2)
            R3.append(r3)
            Eps.append(ep)
            The.append(th)

            x.append(self.consts.c / x0[i-run+1] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                   # [GHz]
            epsilon0, theta0,  = [], []

            for j in range(len(x0[i-run+1])):
                epsilon0.append(R1[i-run+1][j]/ (2 * np.pi * scipy.special.jv(1, 2.405) * R3[i-run+1][j]))                                              # Ellipticity: [millirad]
                theta0.append(R2[i-run+1][j] * 1e3/ (2 * np.pi * scipy.special.jv(2, 2.405) * R3[i-run+1][j] * np.sqrt(1 - 4 * epsilon0[j]**2)))        # Rotation: [millirad]

            # ax.plot(x[i-run+1], epsilon0, label=r'$\epsilon_\text{cell}$' if i-run+1==0 else r'$\epsilon\text{air}$')
            ax.plot(x[i-run+1], theta0, label=r'$\theta_\text{cell}$' if i-run+1==0 else r'$\theta_\text{air}$')

        """
        Subtracting measured background R values from sample R values
        """
        yR1, yR2, yR3, epsilon1, epsilon2, theta1, theta2 = [], [], [], [], [], [], []
        for xval in x0[1]:
            cidx = np.argmin(np.abs(x0[0] - xval))
            # Check if cidx is within the bounds of x0[0]
            if 0 <= cidx < len(x0[0]):
                # Use cidx to access the corresponding elements in y1[0] and y1[1]
                yR1_diff = R1[0][cidx] - R1[1][cidx]
                yR2_diff = R2[0][cidx] - R2[1][cidx]
                yR3_diff = R3[0][cidx] - R3[1][cidx]
                ep_diff = Eps[0][cidx] - Eps[1][cidx]
                th_diff = The[0][cidx] - The[1][cidx]
                yR1.append(yR1_diff)
                yR2.append(yR2_diff)
                yR3.append(yR3_diff)
                epsilon1.append(ep_diff * 1e3)
                theta1.append(th_diff * 1e3)
            # else:
                # # Handle out-of-bounds cases as needed
                # y6.append(np.nan)  # For example, append NaN if cidx is out of bounds
                # y7.append(np.nan)  # For example, append NaN if cidx is out of bounds
                # y8.append(np.nan)  # For example, append NaN if cidx is out of bounds

        x0.append(self.consts.c / x0[1] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                            # [GHz]

        for k in range(len(x0[2])):
            epsilon2.append(yR1[k]/ ( 2 * np.pi * scipy.special.jv(1,2.405) * yR3[k]))                                                                  # Ellipticity: [rad]
            theta2.append(yR2[k] * 1e3 / (2 * np.pi * scipy.special.jv(2,2.405) * yR3[k] * np.sqrt(1 - 4 * epsilon2[k]**2)))                            # Rotation: [rad]
        
        # ax.plot(x0[2], epsilon1, label=r'$\epsilon_\text{cell}$-$\epsilon_\text{air}$')
        # ax.plot(x0[2], epsilon2, label=r'$\frac{\text{R}^\text{cell}_{\omega} - \text{R}^\text{air}_{\omega}}{\text{R}^\text{cell}_\text{dc} - \text{R}^\text{air}_\text{dc}}$')
        ax.plot(x0[2], theta1, label=r'$\theta_\text{cell}$-$\theta_\text{air}$')
        ax.plot(x0[2], theta2, label=r'$\frac{\text{R}^\text{cell}_{2\omega} - \text{R}^\text{air}_{2\omega}}{\text{R}^\text{cell}_\text{dc} - \text{R}^\text{air}_\text{dc}}$')

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        # plt.ylabel(r'Ellipticity (millirad.)', fontsize=25)
        plt.ylabel(r'Faraday Rotation (millirad.)', fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        # plt.title(f'Ellipticity vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.title(f'Faraday Rotation vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} $\mu$W @{date}', fontsize=25)
        # plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity(R)_vs_Frequency_{date}_run{run}-{run+1}.png'))
        plt.savefig(os.path.join(Plots, f'{date}', f'FR(R)_vs_Frequency_{date}_run{run}-{run+1}.png'))
        plt.show()

plotter = Plot()
date_input = '03-21-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.Ellipticity_FR_X(Bristol_path, Lockins_path, 7, 5, -5.103, 860)
# plotter.Ellipticity_FR_R(Bristol_path, Lockins_path, 17, 5, 5.103, 3)
# plotter.theory_plot(0.0718, 5.103, 26)