import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory_calculations import Theory
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
        self.reader = Read()
        self.analyzer = Analyze()

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run+1)

    def plot_process(self, t, X, Y, R, run, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in self.number_of_runs(run):

            if name == '1f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == '2f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'dc':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'm2f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            
            # Apply the scaling factor to X and Y arrays
            X[i] = X[i] * scale_factor
            Y[i] = Y[i] * scale_factor
            R[i] = R[i] * scale_factor


            if i == run-1:
                # Plot X and Y with the appropriate labels
                label_x = (r'$\text{X}_\text{f}$' if name == '1f' else 
                        r'$\text{X}_\text{2f}$' if name == '2f' else 
                        r'$\text{X}_\text{dc}$' if name == 'dc' else
                        r'$\text{X}_\text{m2f}$')
                label_y = (r'$\text{Y}_\text{f}$' if name == '1f' else 
                        r'$\text{Y}_\text{2f}$' if name == '2f' else 
                        r'$\text{Y}_\text{dc}$' if name == 'dc' else
                        r'$\text{Y}_\text{m2f}$')

                label_R = (r'$\text{R}_\text{1f}$' if name == '1f' else 
                        r'$\text{R}_\text{2f}$' if name == '2f' else 
                        r'$\text{R}_\text{dc}$' if name == 'dc' else
                        r'$\text{R}_\text{m2f}$')

                # ax.plot(t[i][8:]/60, X[i][8:], label=label_x, color='r')
                ax.plot(t[i][8:]/60, X[i][8:],  label=label_x, color='r', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=100)
                # ax.plot(t[i][8:]/60, Y[i][8:], label=label_y, color='b')
                ax.plot(t[i][8:]/60, Y[i][8:],  label=label_y, color='b', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=100)
                # ax.scatter(t[i][8:]/60, R[i][8:], label=label_R, color='black', s=50)
                ax.plot(t[i][8:]/60, R[i][8:], label=label_R, color='black', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=100)
            else:
                # ax.plot(t[i][8:]/60, X[i][8:], label=label_x, color='r')
                ax.plot(t[i][8:]/60, X[i][8:],  label=label_x, color='r', linestyle='-', linewidth=1, 
                        marker='x', markersize=10, markevery=100)
                # ax.plot(t[i][8:]/60, Y[i][8:], label=label_y, color='b')
                ax.plot(t[i][8:]/60, Y[i][8:],  label=label_y, color='b', linestyle='-', linewidth=1, 
                        marker='x', markersize=10, markevery=100)
                # ax.scatter(t[i][8:]/60, R[i][8:], label=label_R, color='black', s=50)
                ax.plot(t[i][8:]/60, R[i][8:], label=label_R, color='black', linestyle='-', linewidth=1, 
                        marker='x', markersize=10, markevery=100)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        # plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        # plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(plots, f'{date}', f'{name}_vs_time_{date}_run{run}.png'))
        plt.show()

    def XYR_vs_time(self, lockins_path, name, run, power):
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xm2f, Ym2f = self.reader.read_lockins(lockins_path)
        para, lockins_t, R1f, R2f, Rdc, Rm2f = self.analyzer.R_lockins(lockins_path)
        B_ave, B_variration = self.B_field()

        if name == '1f':
            self.plot_process(lockins_t, X1f, Y1f, R1f, run, name, 'Time (min)',
                                  r'Voltage (mV)', r'1f Voltages vs Time, run' + f'{run}' + 
                                  f', $P$={power} μW') #, $B_z$={B_ave:.3f}$\pm${B_variration:.3f} G, $T$={temp}°C
        elif name == '2f':
            self.plot_process(lockins_t, X2f, Y2f, R2f, run, name, 'Time (min)',
                                  r'Voltage (mV)', r'2f Voltages vs Time, run' + f'{run}' + 
                                  f', $P$={power} μW')
        elif name == 'dc':
            self.plot_process(lockins_t, Xdc, Ydc, Rdc, run, name, 'Time (min)',
                                  r'Voltage (mV)', r'DC Voltages vs Time, run' + f'{run}' + 
                                  f', $P$={power} μW')
        elif name == 'm2f':
            self.plot_process(lockins_t, Xm2f, Ym2f, Rm2f, run, name, 'Time (min)',
                                  r'Voltage (mV)', r'M2f Voltages vs Time, run' + f'{run}' + 
                                  f', $P$={power} μW')

    def B_field(self):
        B_max = np.array([5.2934, 5.2915, 5.2932, 5.2927, 5.2935])
        B_min = np.array([5.2848, 5.2834, 5.2842, 5.2839, 5.2837])

        # Compute the average field
        B_avg = np.round(0.5 * (np.mean(B_max) + np.mean(B_min)),3)
        print("Average Magnetic Field:", B_avg)

        # Compute mean absolute deviation (MAD)
        B_spread = np.round(0.5 * (np.mean(np.abs(B_max - B_avg)) + np.mean(np.abs(B_min - B_avg))),3)
        print("Average Magnetic Field Spread (Variation):", B_spread)

        return B_avg, B_spread
    
if __name__ == "__main__":
    dir_path = os.path.join(
    os.path.expanduser('~'),  # Directory path on personal computer
    'OneDrive', 
    'Files', 
    'Graduate_study', 
    'Research', 
    'PhD_project', 
    'Faraday_rotation_measurements'
    )
    # dir_path = os.path.join(os.getcwd(),  # Directory path on Faraday lab computer
    # 'Faraday_rotation_measurements', 
    # )
    K_vapor = os.path.join(dir_path, 'K_vapor_cell')
    lockins = os.path.join(K_vapor, 'Lockins_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '02-18-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    lockins_path = glob.glob(os.path.join(lockins, date, '*.lvm'))
    plotter.XYR_vs_time(lockins_path, 'm2f', 1, 200)