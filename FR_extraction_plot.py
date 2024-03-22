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

# dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
dir_path = os.path.join(os.getcwd(), 'Faraday Rotation Measurements')
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

    def Ellipticity_vs_Frequency(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        B_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc = self.reader.lockins(lockin_path)
        para, lockins_t, R1f, R2f, Rdc = self.analyzer.R_lockins(lockin_path)
        epsilon, theta = self.analyzer.FR_double_Kvapor(lockin_path)
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        x, x0, X1, X2, X3, R1, R2, R3, Eps, The = [], [], [], [], [], [], [], [], [], []

        for i in range(run-1, run+1):
            B_t[i], Lambda[i] = self.analyzer.filter_data(B_t[i], Lambda[i])

            # trim R values
            B_t[i], Lambda[i], lockins_t[i], X1f[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], X1f[i])
            B_t[i], Lambda[i], lockins_t[i], X2f[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], X2f[i])
            B_t[i], Lambda[i], lockins_t[i], Xdc[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], Xdc[i])
            B_t[i], Lambda[i], lockins_t[i], R1f[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], R1f[i])
            B_t[i], Lambda[i], lockins_t[i], R2f[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], R2f[i])
            B_t[i], Lambda[i], lockins_t[i], Rdc[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], Rdc[i])
            B_t[i], Lambda[i], lockins_t[i], epsilon[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], epsilon[i])
            B_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(B_t[i], Lambda[i], lockins_t[i], theta[i])

            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(B_t[i], lockins_t[i], para[i][2], n)

            # average R values
            Lambd, x1 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], X1f[i][l_idx])
            Lambd, x2 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], X2f[i][l_idx])
            Lambd, x3 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], Xdc[i][l_idx])
            Lambd, r1 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R1f[i][l_idx])
            Lambd, r2 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], R2f[i][l_idx])
            Lambd, r3 = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], Rdc[i][l_idx])
            Lambd, ep = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], epsilon[i][l_idx])
            Lambd, th = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])

            x0.append(Lambd)
            X1.append(x1)
            X2.append(x2)
            X3.append(x3)
            R1.append(r1)
            R2.append(r2)
            R3.append(r3)
            Eps.append(ep)
            The.append(th)

            x.append(self.consts.c / x0[i-run+1] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                   # [GHz]
            epsilon0, theta0,  = [], []

            for j in range(len(x0[i-run+1])):
                epsilon0.append(R1[i-run+1][j] * 1e3/ (2 * np.pi * scipy.special.jv(1, 2.405) * R3[i-run+1][j]))                                             # Ellipticity: [millirad]
                theta0.append(R2[i-run+1][j]/ (2 * np.pi * scipy.special.jv(2, 2.405) * R3[i-run+1][j] * np.sqrt(1 - 4 * epsilon0[j]**2)))       # Rotation: [millirad]

            ax.plot(x[i-run+1], epsilon0, label=r'$\epsilon_7$' if i-run+1==0 else r'$\epsilon_8$')
            # ax.plot(x[i-run+1], theta0, label=r'$\theta_7$' if i-run+1==0 else r'$\theta_8$')

        """
        Subtracting measured background R values from sample R values
        """
        yX1, yX2, yX3, yR1, yR2, yR3, epsilon1, epsilon2, epsilon3, theta1, theta2, theta3 = [], [], [], [], [], [], [], [], [], [], [], []
        for xval in x0[1]:
            cidx = np.argmin(np.abs(x0[0] - xval))
            # Check if cidx is within the bounds of x0[0]
            if 0 <= cidx < len(x0[0]):
                # Use cidx to access the corresponding elements in y1[0] and y1[1]
                yX1_diff = X1[0][cidx] - X1[1][cidx]
                yX2_diff = X2[0][cidx] - X2[1][cidx]
                yX3_diff = X3[0][cidx] - X3[1][cidx]
                yR1_diff = R1[0][cidx] - R1[1][cidx]
                yR2_diff = R2[0][cidx] - R2[1][cidx]
                yR3_diff = R3[0][cidx] - R3[1][cidx]
                ep_diff = Eps[0][cidx] - Eps[1][cidx]
                th_diff = The[0][cidx] - The[1][cidx]
                yX1.append(yX1_diff)
                yX2.append(yX2_diff)
                yX3.append(yX3_diff)
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

        x0.append(self.consts.c / x0[1] * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                                # [GHz]

        for k in range(len(x0[2])):
            print(yX1[k])
            epsilon2.append(yX1[k] * 1e3 / ( 2 * np.pi * scipy.special.jv(1,2.405) * yX3[k]))                                                       # Ellipticity: [rad]
            epsilon3.append(yR1[k] * 1e3 / ( 2 * np.pi * scipy.special.jv(1,2.405) * yR3[k]))                                                       # Ellipticity: [rad]
            theta2.append(yX2[k] * 1e3 / (2 * np.pi * scipy.special.jv(2,2.405) * yX3[k] * np.sqrt(1 - 4 * epsilon2[k]**2)))                   # Rotation: [rad]
            theta3.append(yR2[k] * 1e3 / (2 * np.pi * scipy.special.jv(2,2.405) * yR3[k] * np.sqrt(1 - 4 * epsilon2[k]**2)))                   # Rotation: [rad]
        
        ax.plot(x0[2], epsilon1, label=r'$\epsilon_7$-$\epsilon_8$')
        ax.plot(x0[2], epsilon2, label=r'$X_7$-$X_8$')
        ax.plot(x0[2], epsilon3, label=r'$R_7$-$R_8$')
        
        # ax.plot(x0[2], theta1, label=r'$\theta_7$-$\theta_8$')
        # ax.plot(x0[2], theta2, label=r'$X_7$-$X_8$')
        # ax.plot(x0[2], theta3, label=r'$R_7$-$R_8$')
        

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.ylabel(r'Ellipticity (millirad.)', fontsize=25)
        # plt.ylabel(r'Faraday Rotation (millirad.)', fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(f'Ellipticity vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        # plt.title(f'Raraday Rotation vs Frequency, run{run}-{run+1}, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'Ellipticity_vs_Frequency_{date}_run{run}-{run+1}.png'))
        # plt.savefig(os.path.join(Plots, f'{date}', f'FR_vs_Frequency_{date}_run{run}-{run+1}.png'))
        plt.show()

    def FR_vs_Frequency(self, lambda_path, lockin_path, run, n, B, power):
        """
        Plot function of measured FR angle vs Frequency
        """
        Bristol_t, Lambda = self.reader.Bristol(lambda_path)
        para, lockins_t, R1f, R2f, Rdc = self.analyzer.R_lockins(lockin_path)
        epsilon, theta = self.analyzer.FR_double_Kvapor(lockin_path)
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        x0, y0 = [], []
        for i in range(run-1, run+1):
            Bristol_t[i], Lambda[i] = self.analyzer.filter_data(Bristol_t[i], Lambda[i])
            Bristol_t[i], Lambda[i], lockins_t[i], theta[i] = self.analyzer.trim_data(Bristol_t[i], Lambda[i], lockins_t[i], theta[i])
            l_idx, b_idx = self.analyzer.calculate_interval_and_indices(Bristol_t[i], lockins_t[i], para[i][2], n)
            Lambd, Thet = self.analyzer.calculate_averages(b_idx, Lambda[i], Lambda[i][b_idx], theta[i][l_idx])
            x0.append(self.consts.c / Lambd * 1e-9 - self.consts.Nu39_D2 * 1e-9)                                                   # [GHz]
            y0.append(Thet[1:] * 1e3)                                                                                              # [millirad]
        print(len(x0[0]), len(x0[1]))
        print(len(y0[0]), len(y0[1]))
        y1 = []
        for xval in x0[1]:
            cidx = np.argmin(np.abs(x0[0] - xval))
            y1.append(y0[0][cidx] - y0[1][cidx] )
        ax.plot(x0[0], y1, label=r'L->H Wide Scan')

        # y1 = []
        # for xval in x0[1]:
        #     cidx = np.argmin(np.abs(x0[3] - xval))
        #     y1.append(y0[3][cidx] - y0[1][cidx] )
        
        # ax.plot(x0[1], y1, label=r'H->L Wide Scan')
        
        # lambda_theo = np.linspace(766.69*1e-9, 766.71*1e-9, 2000)
        # y_theo = self.theory.FR_theta1(lambda_theo, 0.0718, B*1e-4, 26, self.consts.Lambda39_D1, self.consts.Lambda39_D2)
        # x_theo = self.consts.c / lambda_theo * 1e-9 - self.consts.Nu39_D2  * 1e-9
        # plt.plot(x_theo, y_theo * 1e3, '--', color='red', label='Theory')

        plt.xlabel(r'Frequency (GHz)', fontsize=25)
        plt.ylabel(r'Faraday Rotation (millirad.)', fontsize=25)
        plt.xticks(np.arange(-5, 7, 1), fontsize=25)
        plt.yticks(fontsize=25)
        # plt.xlim(-4,6)
        # plt.ylim(-0.15,0.25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(f'Faraday Rotation vs Frequency, $B_z$={B} G, $P$={power} nW @{date}', fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'FR_vs_Frequency_{date}_run{run}-{run+3}.png'))
        plt.show()

    def theory_plot(self, l, B, T):
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        x1 = np.linspace(760*1e-9, 777*1e-9, 20000)
        y1 = self.theory.FR_theta1(x1, l, B*1e-4, T, self.consts.Lambda39_D1, self.consts.Lambda39_D2)
        # plt.plot(x1 * 1e9, y1 * 1e3, '-', color='red', label='Theory1')

        x2 = np.linspace(self.consts.c/(777*1e-9), self.consts.c/(760*1e-9), 20000)
        y2 = self.theory.FR_theta2(x2, l, B*1e-4, T, self.consts.Nu39_D1, self.consts.Nu39_D2)
        plt.plot(x2 * 1e-9, y2 * 1e3, '--', color='blue', label='Theory1')
        plt.xticks(np.arange(min(x2*1e-9), max(x2*1e-9), (max(x2*1e-9)-min(x2*1e-9))/10))
        ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        # plt.ylim(0, 2e5)
        plt.grid(True)
        plt.show()

plotter = Plot()
date_input = '03-21-2024'
date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
plotter.Ellipticity_vs_Frequency(Bristol_path, Lockins_path, 7, 5, 5.103, 860)
# plotter.FR_vs_Frequency(Bristol_path, Lockins_path, 7, 5, 5.103, 860)
# plotter.theory_plot(0.0718, 5.103, 26)