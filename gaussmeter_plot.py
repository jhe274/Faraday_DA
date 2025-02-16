import os, glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from data_reader import DataReader as Read
from data_analyzer import DataAnalyzer as Analyze

class Plot:

    def __init__(self):
        self.reader = Read()
        self.analyzer = Analyze()

    def number_of_runs(self, run):
        """
        Determine the range of runs to analyze.
        :param run: Current run index
        :return: A range of run indices
        """
        return range(run-1, run)
    
    def process_temperature(self, temps, run):
        T_mean = np.mean(temps[run-1])
        T_std = np.std(temps[run-1], ddof=1) / np.sqrt(len(temps[run-1]))

        return T_mean, T_std
    
    def plot_process(self, x, y1, y2, run, xlabel, ylabel_y1, ylabel_y2, file_name):
        # Create a figure and axes for plotting
        fig, ax1 = plt.subplots(1, 1, figsize=(25.60, 14.40))
        # Create second Y-axis
        ax2 = ax1.twinx()  # Create a second y-axis that shares the same x-axis

        for i in self.number_of_runs(run):
            if i == run-1:
                y1_fit, slope, intercept, y1_mean, residuals, y1_std = self.analyzer.drift_fit(x[i]/60, y1[i])
                T_mean, T_std = self.process_temperature(y2, run)
                
                lower_bound = slope*x[i]/60 + intercept - y1_std
                upper_bound = slope*x[i]/60 + intercept + y1_std

                ax1.plot(x[i]/60, y1[i], color='C0', alpha=0.4, label=f'$\\overline{{B_z}}$={round(y1_mean,3):.3f} G')
                ax1.plot(x[i]/60, y1_fit, '--', color='b', label=fr'$\nabla_t B_z$={round(slope*1e3,1):.1f} mG/min')
                ax1.fill_between(x[i]/60, lower_bound, upper_bound, facecolor='C0', alpha=0.4, label=f'$\\sigma$={round(y1_std*1e3)} mG')
                # ax1.fill_between(x[i]/60, upper_bound, y1[i], where=y1[i] > upper_bound, fc='red', alpha=0.4, interpolate=True)
                # ax1.fill_between(x[i]/60, lower_bound, y1[i], where=y1[i] < lower_bound, fc='red', alpha=0.4, interpolate=True)

                ax1.set_xlabel(xlabel, fontsize=25)
                ax1.set_ylabel(ylabel_y1, fontsize=25)
                ax1.tick_params(axis='x', labelsize=25)
                ax1.tick_params(axis='y', labelsize=25)

                ax2.plot(x[i]/60, y2[i], label=f'$\\overline{{T}}$={round(T_mean,2):.2f}±{round(T_std,2):.2f}°C', color='black', linestyle='-', linewidth=1, 
                        marker='^', markersize=10, markevery=500)
                ax2.set_ylabel(ylabel_y2, fontsize=25)
                ax2.tick_params(axis='y', labelsize=25)
            else:
                y1_fit, slope, intercept, y1_mean, residuals, y1_std = self.analyzer.drift_fit(x[i]/60, y1[i])
                
                # ax1.plot(x[i]/60, y1[i], label=r'$B_0$, run2', color='r', alpha=0.4, linestyle='-', linewidth=1, 
                #         marker='x', markersize=10, markevery=200)
                # ax1.plot(x[i]/60, y1_fit, label=f'Linear regression: run{run+1}', color='r', linestyle='--', linewidth=1)

                lower_bound = slope*x[i]/60 + intercept - y1_std
                upper_bound = slope*x[i]/60 + intercept + y1_std

                ax1.plot(x[i]/60, y1[i], color='C3', alpha=0.4, label=f'$\\overline{{B_z}}$={round(y1_mean,3):.3f} G')
                ax1.plot(x[i]/60, y1_fit, '--', color='r', label=fr'$\nabla_t B_z$={round(slope*1e3,3):.3f} mG/min')
                ax1.fill_between(x[i]/60, lower_bound, upper_bound, facecolor='C3', alpha=0.4, label=f'$\\sigma$={round(y1_std*1e3)} mG')
                # ax1.fill_between(x[i]/60, upper_bound, y1[i], where=y1[i] > upper_bound, fc='r', alpha=0.4, interpolate=True)
                # ax1.fill_between(x[i]/60, lower_bound, y1[i], where=y1[i] < lower_bound, fc='r', alpha=0.4, interpolate=True)

                ax2.plot(x[i]/60, y2[i], label=f'$\\overline{{T}}$={round(T_mean,2):.2f}±{round(T_std,2):.2f}°C', color='black', linestyle='-', linewidth=1, 
                        marker='x', markersize=10, markevery=500)

            # plt.title(f'run{run}', fontsize=25)
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=20)
        plt.grid(False)
        save_path = os.path.join(plots, date, file_name)
        plt.savefig(save_path)
        plt.show()
    
    def gaussmter_vs_time(self, gaussmeter_path, run):
        timestamps, B0s, temps = self.reader.read_gaussmeter(gaussmeter_path)

        self.plot_process(timestamps, B0s, temps, run, 'Time (min)',
                            r'Magnetic flux density (G)', r'Temperature (°C)', 
                            f'Magnetic_field_and_temperature_{date}_run{run}.png')
    

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
    gaussmeter = os.path.join(K_vapor, 'Gaussmeter_data')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '02-13-2025'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    gaussmeter_path = glob.glob(os.path.join(gaussmeter, date, '*.csv'))
    plotter.gaussmter_vs_time(gaussmeter_path, 7)