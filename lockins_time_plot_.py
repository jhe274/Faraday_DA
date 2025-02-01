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
        return range(run-1, run)

    def XYplot(self, t, X, Y, run, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))
        for i in self.number_of_runs(run):
            if name == 'f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == '2f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'dc':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'mod':
                scale_factor = 1e3  # RMS Voltage: [mV]
            
            # Apply the scaling factor to X and Y arrays
            X[i] = X[i] * scale_factor
            Y[i] = Y[i] * scale_factor

            # Plot X and Y with the appropriate labels
            label_x = (r'$\text{X}_\text{f}$' if name == 'f' else 
                    r'$\text{X}_\text{2f}$' if name == '2f' else 
                    r'$\text{X}_\text{dc}$' if name == 'dc' else
                    r'$\text{Y}_\text{mod}$')
            label_y = (r'$\text{Y}_\text{f}$' if name == 'f' else 
                    r'$\text{Y}_\text{2f}$' if name == '2f' else 
                    r'$\text{Y}_\text{dc}$' if name == 'dc' else
                    r'$\text{Y}_\text{mod}$')

            ax.plot(t[i], X[i], color='r', label=label_x, linestyle='-', linewidth=1, marker='^', markevery=5, markersize=10)
            ax.plot(t[i], Y[i], color='b', label=label_y, linestyle='-', linewidth=1, marker='x', markevery=5, markersize=10)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        # plt.xticks(np.arange(-5, 6, 1), fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        # plt.grid(True)
        ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'XY{name}_{date}_run{i}-{i+1}.pdf'))
        plt.show()

    def Rplot(self, t, R, run, name, xlabel, ylabel, title):
        fig, ax = plt.subplots(1, 1, figsize=(25, 12))
        for i in self.number_of_runs(run):
            if name == 'f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == '2f':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'dc':
                scale_factor = 1e3  # RMS Voltage: [mV]
            elif name == 'mod':
                scale_factor = 1e3  # RMS Voltage: [mV]
            
            # Apply the scaling factor to X and Y arrays
            R[i] = R[i] * scale_factor

            label_R = (r'$\text{R}_\text{f}$' if name == 'f' else 
                    r'$\text{R}_\text{2f}$' if name == '2f' else 
                    r'$\text{R}_\text{dc}$' if name == 'dc' else
                    r'$\text{R}_\text{mod}$')
            ax.scatter(t[i], R[i], label=label_R, color='r', s=10)

        plt.xlabel(xlabel, fontsize=25)
        plt.ylabel(ylabel, fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        # ax.get_xaxis().set_major_formatter(plt.FormatStrFormatter('%.3f'))
        # plt.grid(True)
        # ax.legend(loc='best', fontsize=25)
        plt.title(title, fontsize=25)
        plt.savefig(os.path.join(Plots, f'{date}', f'R{name}_{date}_run{i}-{i+1}.png'))
        # plt.savefig(os.path.join(Plots, f'{date}', f'{name}_{date}_run{i}.png'))
        plt.show()

    def XY_vs_time(self, lockins_path, name, run, B, power):
        para, lockins_t, X1f, Y1f, X2f, Y2f, Xdc, Ydc, Xmod, Ymod = self.reader.read_lockins(lockins_path)
        if name == 'f':
            self.XYplot(lockins_t, X1f, Y1f, run-1, name, 'Time (s)',
                                  r'$\text{XY}_{1f}$ (mV)', r'$\text{XY}_{1f}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
        elif name == '2f':
            self.XYplot(lockins_t, X2f, Y2f, run-1, name, 'Time (s)',
                                  r'$\text{XY}_{2f}$ (mV)', r'$\text{XY}_{2f}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
        elif name == 'dc':
            self.XYplot(lockins_t, Xdc, Ydc, run-1, name, 'Time (s)',
                                  r'$\text{XY}_\text{dc}$ (mV)', r'$\text{XY}_\text{dc}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
        elif name == 'mod':
            self.XYplot(lockins_t, Xmod, Ymod, run-1, name, 'Time (s)',
                                  r'$\text{XY}_\text{mod}$ (mV)', r'$\text{XY}_\text{mod}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
            
    def R_vs_time(self, lockins_path, name, run, B, power):
        para, lockins_t, R1f, R2f, Rdc, Rmod = self.analyzer.R_lockins(lockins_path)

        if name == 'f':
            self.Rplot(lockins_t, R1f, run-1, name, 'Time (s)',
                                  r'$\text{R}_{1f}$ (mV)', r'$\text{R}_{1f}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
        elif name == '2f':
            self.Rplot(lockins_t, R2f, run-1, name, 'Time (s)',
                                  r'$\text{R}_{2f}$ (mV)', r'$\text{R}_{2f}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
        elif name == 'dc':
            self.Rplot(lockins_t, Rdc, run-1, name, 'Time (s)',
                                  r'$\text{R}_\text{dc}$ (mV)', r'$\text{R}_\text{dc}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))
        elif name == 'mod':
            self.Rplot(lockins_t, Rmod, run-1, name, 'Time (s)',
                                  r'$\text{R}_\text{mod}$ (mV)', r'$\text{R}_\text{mod}$ vs Time, run' + f'{run}-{run+1}' + 
                                  f', $B_z$={B} G, $P$={power} μW' + ' @'+ str(date))

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
    Bristol = os.path.join(K_vapor, 'Bristol_data')
    Lockins = os.path.join(K_vapor, 'Lockins_data')
    Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '01-31-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
    Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
    # plotter.XY_vs_time(Lockins_path, 'mod', 1, (5.283,5.293), 270)
    plotter.R_vs_time(Lockins_path, 'mod', 1, (5.283,5.293), 270) 