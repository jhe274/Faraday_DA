import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Sample data (replace this with your actual data)
wavelength = np.array([400, 450, 500, 550, 600, 650, 700])
voltage = np.array([1.2, 1.5, 2.0, 2.3, 2.5, 2.8, 3.0])

# Define the function for the baseline fit (you can adjust the degree)
def baseline_fit(x, *params):
    return np.polyval(params, x)

# Initial guess for the parameters (based on your data)
initial_guess = [1.0, 0.0]  # Adjust the values based on your data and degree of the polynomial

# Perform the curve fit
params, covariance = curve_fit(baseline_fit, wavelength, voltage, p0=initial_guess)

# Generate the fitted baseline
fit_baseline = baseline_fit(wavelength, *params)

# Plot the original data and the fitted baseline
plt.scatter(wavelength, voltage, label='Original Data')
plt.plot(wavelength, fit_baseline, color='red', label='Baseline Fit')

plt.xlabel('Wavelength')
plt.ylabel('Photodiode Output Voltage')
plt.legend()
plt.show()

# Laser power without atomic absorption is related to the baseline fit
laser_power_without_absorption = baseline_fit(some_wavelength_value, *params)
print("Estimated Laser Power without Atomic Absorption:", laser_power_without_absorption)
