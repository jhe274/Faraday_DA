import os, glob
import datetime as dt
from constants import Constants as Consts
from theory import Theory
from read import Read
from write import Write
from analyze import Analyze
from plot import Plot
import matplotlib.pyplot as plt

cwd = os.getcwd()
dir_path = os.path.join(cwd, 'Research', 'PhD Project', 'Faraday Rotation Measurements', 'K vapor cell')
Bristol = os.path.join(dir_path, 'Bristol data')
TC300 = os.path.join(dir_path, 'TC300 data')
TriMod = os.path.join(dir_path, 'Triple modulation data')
Analysis = os.path.join(dir_path, 'Data Analysis')



def write_theta_vs_lambda(i, n):
    Write.Theta_and_lambda_data_to_txt(date, Bristol_path, TriMod_path, i, n)

def plot_theta_theory_vs_lambda():
    # Plot.theta_theory_vs_lambda(1.91*pow(10,-2), 12000*pow(10,-4))
    Plot.theta_theory_vs_lambda((7.5-0.159*2)*1e-2, -9.68*1e-4, range(40,50,10), Consts.Lambda_D1, Consts.Lambda_D2)

def plot_temp_vs_time(i):
    Plot.temp_vs_time(TC300_path, TriMod_path, i)

def plot_lambda_vs_time(i, n):
    Plot.lambda_vs_time(Bristol_path, TriMod_path, i, n)

def plot_Vdc_vs_time(i):
    Plot.Vdc_vs_time(TriMod_path, i)

def plot_V2f_vs_time(i):
    Plot.V2f_vs_time(TriMod_path, i)

def plot_Vmod_vs_time(i):
    Plot.Vmod_vs_time(TriMod_path, i)

def plot_theta_mea_vs_time(i, n):
    Plot.theta_mea_vs_time(Bristol_path, TriMod_path, i, n)

def plot_theta_mea_vs_lambda(i, n):
    Plot.theta_mea_vs_lambda(Bristol_path, TriMod_path, i, n)

def plot_theta_mea_vs_nu(i, n):
    Plot.theta_mea_vs_nu(Bristol_path, TriMod_path, i, n)
    plt.savefig(os.path.join(Analysis , 'Plots', f'{date}', f'Theta_vs_nu_{date}_run{i+1}.png'))

def plot_combined_theta_vs_lambda():
    Plot.combined_theta_vs_lambda()

# Function to get date and file paths
def get_date_and_paths():
    date_input = input("Enter the date (MM-DD-YYYY): ")
    date = dt.datetime.strptime(date_input, '%m-%d-%Y').strftime('%m-%d-%Y')
    Bristol_path = glob.glob(os.path.join(Bristol, date, '*.csv'))
    TC300_path = glob.glob(os.path.join(TC300, date, '*.csv'))
    TriMod_path = glob.glob(os.path.join(TriMod, date, '*.lvm'))
    return date, Bristol_path, TC300_path, TriMod_path

# Menu options
options = {
    "1": write_theta_vs_lambda,
    "2": plot_theta_theory_vs_lambda,
    "3": plot_temp_vs_time,
    "4": plot_lambda_vs_time,
    "5": plot_Vdc_vs_time,
    "6": plot_V2f_vs_time,
    "7": plot_Vmod_vs_time,
    "8": plot_theta_mea_vs_time,
    "9": plot_theta_mea_vs_lambda,
    "10": plot_theta_mea_vs_nu,
    "11": plot_combined_theta_vs_lambda
}

# Main menu
print("---------------------------------------------------------------------------------------------------")
print("Type '1' for writting measured polairzation rotation, wavelength to csv file")
print("Type '2' for plotting theoretical polarization rotation vs wavelength")
print("Type '3' for plotting measured temperatures of potassium vapor cell over time")
print("Type '4' for plotting measured wavelength over time")
print("Type '5' for plotting measured Vdc over time")
print("Type '6' for plotting measured V2f over time")
print("Type '7' for plotting measured Vmod over time")
print("Type '8' for plotting measured polarization rotation over time")
print("Type '9' for plotting measured polarization rotation vs wavelength")
print("Type '10' for plotting measured polarization rotation vs frequency detuning")
print("Type '11' for plotting combined polarization rotations vs wavelengths")
print("---------------------------------------------------------------------------------------------------")
option = input("Option: ")

if option in options:
    try:
        selected_function = options[option]
        if selected_function == plot_theta_theory_vs_lambda:
            plot_theta_theory_vs_lambda()
        elif selected_function == plot_combined_theta_vs_lambda:
                plot_combined_theta_vs_lambda()
        else:
            date, Bristol_path, TC300_path, TriMod_path = get_date_and_paths()
            i = input("Select the number of run: ")
            i = int(i) - 1
            if selected_function == plot_temp_vs_time:
                plot_temp_vs_time(i)
            elif selected_function == plot_Vdc_vs_time:
                plot_Vdc_vs_time(i)
            elif selected_function == plot_V2f_vs_time:
                plot_V2f_vs_time(i)
            elif selected_function == plot_Vmod_vs_time:
                plot_Vmod_vs_time(i)
            else:
                n = input("Select the number of integer multiples of time constants to skip in the measurements: ")
                n = int(n)
                if selected_function == write_theta_vs_lambda:
                    write_theta_vs_lambda(i, n)
                else:
                    selected_function(i, n)
    except ValueError:
        print("Exiting...")
else:
    print("Invalid Option. Exiting...")