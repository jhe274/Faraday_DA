import numpy as np
import matplotlib.pyplot as plt
import allantools

# Set the parameters for synthetic data
num_points = 173  # Number of data points
T = 3600  # Total time in seconds (1 hour)
time = np.linspace(0, T, num_points)  # Time axis

# Generate synthetic white noise (random Gaussian fluctuations)
white_noise_std = 0.6e-6  # Standard deviation of 0.6 μrad
white_noise = np.random.normal(0, white_noise_std, len(time))

# Generate slow temperature drift (linear drift)
drift_slope = np.random.uniform(-0.5e-6, 0.5e-6)  # Random drift slope in μrad/h
temperature_drift = drift_slope * (time / 3600)  # Convert seconds to hours

# Total synthetic Faraday rotation data
faraday_rotation = white_noise + temperature_drift

# Compute Allan deviation
taus, adev, _, _ = allantools.oadev(faraday_rotation, rate=num_points/T, data_type="freq")

# Plot the simulated Faraday rotation data
plt.figure(figsize=(10, 4))
plt.plot(time / 3600, faraday_rotation * 1e6, label="Simulated Data", color="b")
plt.xlabel("Time (hours)")
plt.ylabel("Faraday Rotation (μrad)")
plt.title("Synthetic Faraday Rotation Measurement with White Noise and Drift")
plt.legend()
plt.grid()
plt.show()

# Plot the Allan deviation
plt.figure(figsize=(8, 6))
plt.loglog(taus, adev * 1e6, marker="o", linestyle="-", color="r")
plt.xlabel("Averaging Time, τ (seconds)")
plt.ylabel("Allan Deviation (μrad)")
plt.title("Allan Deviation Analysis of Faraday Rotation Data")
plt.grid(which="both", linestyle="--")
plt.show()

