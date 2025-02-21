import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the sine wave function with unknown parameters
def sine_wave(t, offset, amplitude, phase):
    return offset + amplitude * np.sin(2 * np.pi * 0.5 * t + phase)

# Generate synthesized data with some noise
np.random.seed(42)  # For reproducibility
t = np.linspace(0, 10, 100)  # Time vector from 0 to 10 seconds
true_offset = 1.5
true_amplitude = 2.0
true_phase = np.pi / 4  # 45 degrees in radians

# Generate the sine wave with some added noise
y = sine_wave(t, true_offset, true_amplitude, true_phase) + 0.5 * np.random.normal(size=len(t))

# Fit the sine wave to the data
popt, pcov = curve_fit(sine_wave, t, y, p0=[1, 1, 0])

# Extract the fitted parameters
fitted_offset, fitted_amplitude, fitted_phase = popt

# Generate the fitted curve
y_fit = sine_wave(t, *popt)

# Plot the synthesized data and the fitted curve
plt.figure(figsize=(10, 6))
plt.scatter(t, y, label='Synthesized Data with Noise', color='blue')
plt.plot(t, y_fit, label='Fitted Sine Wave', color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Fitting a Sine Wave to Synthesized Data')
plt.legend()
plt.grid(True)
plt.show()

# Print the true and fitted parameters
print(f"True Parameters: Offset={true_offset}, Amplitude={true_amplitude}, Phase={true_phase}")
print(f"Fitted Parameters: Offset={fitted_offset}, Amplitude={fitted_amplitude}, Phase={fitted_phase}")