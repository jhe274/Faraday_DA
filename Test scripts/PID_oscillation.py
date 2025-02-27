import numpy as np
import matplotlib.pyplot as plt

# Simulating PID overtuning effects on laser frequency stability
t = np.linspace(0, 10, 1000)  # Time in seconds

# Normal response (well-tuned PID)
stable_frequency = 0.5 * np.exp(-0.5 * t) * np.cos(2 * np.pi * 0.5 * t)

# Overshoot due to high P gain (P too high, oscillations before settling)
high_P_oscillations = np.exp(-0.1 * t) * np.cos(2 * np.pi * 3 * t)  # Higher frequency oscillations

# Integral windup (I too high, slow recovery with sustained oscillations)
high_I_windup = 0.7 * np.exp(-0.02 * t) * np.cos(2 * np.pi * 1.5 * t)

# Excessive D gain (D too high, amplifies high-frequency noise)
high_D_noise = 0.5 * np.exp(-0.5 * t) * np.cos(2 * np.pi * 0.5 * t) + 0.05 * np.random.randn(len(t))

# Plot the responses
plt.figure(figsize=(10, 6))
plt.plot(t, stable_frequency, label="Well-tuned PID (Stable Lock)", linewidth=2, color="green")
plt.plot(t, high_P_oscillations, label="High P Gain (Oscillatory Lock)", linestyle="dashed", color="red")
plt.plot(t, high_I_windup, label="High I Gain (Windup & Slow Recovery)", linestyle="dotted", color="blue")
plt.plot(t, high_D_noise, label="High D Gain (Amplifies Noise)", linestyle="dashdot", color="purple")

# Labels and legend
plt.xlabel("Time (s)")
plt.ylabel("Laser Frequency Deviation (arbitrary units)")
plt.title("Effect of Overtuned PID Parameters on Laser Frequency Stability")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
