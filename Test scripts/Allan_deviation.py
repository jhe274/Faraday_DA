import numpy as np
import allantools
import matplotlib.pyplot as plt

# Generate synthetic data
np.random.seed(42)  # For reproducibility
measurement_time = 4 * 3600  # 4 hours in seconds
sampling_rate = 1  # 1 sample per second
num_samples = measurement_time // sampling_rate
modulation_frequency = 0.5  # Hz
modulation_amplitude = 1e-6  # Small unknown amplitude in radians

# Time array
time = np.arange(0, measurement_time, sampling_rate)

# Simulated signal: Small sinusoidal modulation with noise
signal = modulation_amplitude * np.sin(2 * np.pi * modulation_frequency * time)
noise = np.random.normal(0, modulation_amplitude / 10, num_samples)  # Add some noise
measured_signal = signal + noise

# Allan deviation analysis
taus, adev, _, _ = allantools.oadev(measured_signal, rate=sampling_rate, data_type='freq')

# Plot Allan deviation
plt.figure(figsize=(8, 6))
plt.loglog(taus, adev, marker='o', linestyle='-', label='Allan Deviation')
plt.xlabel('Averaging Time (s)')
plt.ylabel('Allan Deviation')
plt.title('Allan Deviation Analysis')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()

# Print key results
min_adev = np.min(adev)
optimal_tau = taus[np.argmin(adev)]
print(f"Minimum Allan deviation: {min_adev:.2e} at averaging time {optimal_tau:.1f} s")
