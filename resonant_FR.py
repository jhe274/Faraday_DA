import os, glob
import numpy as np
import matplotlib.pyplot as plt
from constants import Constants as Consts
from theory_calculations import Theory

class Plot:
    """
    A class to handle data visualization, analysis, and processing for Faraday rotation measurements.
    """
    def __init__(self):
        """
        Initialize the Plot class with constants, theory, reader, and analyzer modules.
        """
        self.consts = Consts()  # Load physical constants
        self.theory = Theory()  # Load theoretical models
        self.l = (7.5 - 0.159 * 2) * 1e-2  # Optical path length in meters

    def plot(self, l):
        """
        Plot the Faraday rotation data.
        """
        freqs = np.linspace(self.consts.K39_D2_Hz - 5e9, self.consts.K39_D2_Hz + 5e9, 4000) # [Hz]
        detuning = (freqs - self.consts.K39_D2_Hz) * 1e-9 # [GHz]
        dia_FR = self.theory.diamagnetic_FR(freqs, l, 2, 294, -6, 100, 1000) * 1e6 # [µrad]
        para_FR = self.theory.paramagnetic_FR(freqs, l, 2, 294, 0.01, 100, 1000) * 1e6 # [µrad]

        fig, ax = plt.subplots(1, 1, figsize=(25.60, 14.40))

        ax.plot(detuning, dia_FR, c='r', label='Diamagnetic FR')
        ax.plot(detuning, para_FR, c='b',label='Paramagnetic FR')
        ax.plot(detuning, dia_FR + para_FR, c='black', label='Total FR')

        ax.set_xlabel(r'Frequency Detuning, $\nu$ (GHz)', fontsize=30)
        ax.set_ylabel(r'Faraday Rotation, $\theta$ (µrad)', fontsize=30)
        ax.set_xticks(np.arange(-5, 6, 1))
        ax.tick_params(axis='both', which='major', labelsize=30)
        ax.legend(fontsize=30, loc="best")  # Use custom handler

        plt.tight_layout()
        save_path = os.path.join(plots, 'Kadleck_FR_theory_1_percent_polarization.png')
        plt.savefig(save_path)
        plt.show()

if __name__ == "__main__":
    dir_path = os.path.join(
    os.path.expanduser('~'),  # Expands to your home directory
    'OneDrive', 
    'Files', 
    'Graduate_study', 
    'Research', 
    'PhD_project', 
    'Faraday_rotation_measurements'
    )
    K_vapor = os.path.join(dir_path, 'K_vapor_cell')
    plots = os.path.join(dir_path, 'Data_analysis', 'Plots')
    processed_path = os.path.join(dir_path, 'Data_analysis', 'Processed_data')
    
    plotter = Plot()
    plotter.plot(plotter.l)