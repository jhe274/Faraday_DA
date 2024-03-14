import os
import glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from constants import Constants as Consts
from theory import Theory
from read import Read
from analyze import Analyze

class Plot:

    def __init__(self):
        self.analyzer = Analyze()

    def rdc_vs_rtd(self, lockins_path, B, power, date):
        para, lockins_t, R1f, R2f, Rdc, epsilon, theta = self.analyzer.Double_modu_theta(lockins_path)

        fig, ax = plt.subplots(figsize=(25, 12))
        rtd = np.linspace(2.2, 3, 17)
        x = np.arange(1.8, 3.15, 0.05)

        def process_data(Rdc_slice, label):
            Vdc = [np.average(Rdc_slice[i::len(rtd)]) for i in range(len(rtd))]
            y = [v * 1e3 for v in Vdc]
            coefficients = np.polyfit(rtd, y, 1)
            m, c = coefficients
            fitted_y = m * x + c

            # Plot the line and get the intersection points with previous lines
            ax.plot(x, fitted_y, label=label)
            if ax.lines:
                for line in ax.lines[:-1]:
                    x_intersect, y_intersect = self.find_intersection(x, fitted_y, line.get_xdata(), line.get_ydata())
                    if x_intersect is not None:
                        ax.scatter(x_intersect, y_intersect, color='black', s=100)
                        ax.annotate(f'({x_intersect:.2f}, {y_intersect:.2f})', (x_intersect, y_intersect), textcoords="offset points", xytext=(-20,20), ha='center')

            return y, fitted_y

        labels = [r'$\theta_a$=30°', r'$\theta_a$=75°', r'$\theta_a$=120°']
        colors = ['b', 'g', 'r']
        data_slices = [Rdc[i*len(rtd):(i+1)*len(rtd)] for i in range(len(labels))]

        for label, color, Rdc_slice in zip(labels, colors, data_slices):
            y, fitted_y = process_data(Rdc_slice, label)
            ax.scatter(rtd, y, color=color)
            ax.plot(x, fitted_y, color=color)

        ax.set_xlabel('RTD (rad)', fontsize=25)
        ax.set_ylabel(r'$R_{dc}$ (mV)', fontsize=25)
        plt.xticks(np.arange(1.8, 3.2, 0.1), fontsize=25)
        plt.yticks(fontsize=25)
        ax.grid(True)
        ax.legend(loc='best', fontsize=25)
        ax.set_title(fr'$R_{{dc}}$ vs RTD, $B_z$={B} G, $P$={power} $\mu$W @{date}', fontsize=25)
        plt.savefig(os.path.join(Plots, date, f'PEM_rtd_cali_{date}.png'))
        plt.show()

    def find_intersection(self, x1, y1, x2, y2):
        # Find intersection of two lines given their equations
        coefficients1 = np.polyfit(x1, y1, 1)
        coefficients2 = np.polyfit(x2, y2, 1)

        m1, c1 = coefficients1
        m2, c2 = coefficients2

        # If the lines are parallel, there is no intersection
        if m1 == m2:
            return None, None

        x_intersect = (c2 - c1) / (m1 - m2)
        y_intersect = m1 * x_intersect + c1

        return x_intersect, y_intersect

if __name__ == "__main__":
    dir_path = os.path.join(os.getcwd(), 'Research', 'PhD Project', 'Faraday Rotation Measurements')
    K_vapor = os.path.join(dir_path, 'K vapor cell')
    Bristol = os.path.join(K_vapor, 'Bristol data')
    Lockins = os.path.join(K_vapor, 'Lockins data')
    Plots = os.path.join(dir_path, 'Data_analysis', 'Plots')

    plotter = Plot()
    date_input = '03-13-2024'
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
    Lockins_path = glob.glob(os.path.join(Lockins, date, '*.lvm'))
    plotter.rdc_vs_rtd(Lockins_path, 5.103, 1.5, date)
