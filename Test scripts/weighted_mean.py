import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Seed for reproducibility
np.random.seed(42)

# Synthesized time-series data
time = np.arange(0, 240, 1)  # 4 hours of data in minutes
# Linear drift component
drift = 0.05 * time
# Random noise component
noise = np.random.normal(0, 1, len(time))
# Combined signal
measurements = drift + noise

# Create a DataFrame
df = pd.DataFrame({'Time': time, 'Measurements': measurements})

# Perform linear regression to model the drift
coefficients = np.polyfit(df['Time'], df['Measurements'], deg=1)
trend = np.polyval(coefficients, df['Time'])

# Detrend the data
df['Detrended'] = df['Measurements'] - trend

# Calculate rolling variance
window_size = 30  # 30-minute window
df['Rolling_Var'] = df['Detrended'].rolling(window=window_size, center=True).var()

# To avoid division by zero or NaN values, replace NaNs with a large number
df['Rolling_Var'].fillna(df['Rolling_Var'].max(), inplace=True)
print(df['Rolling_Var'])

# Calculate weights as the inverse of rolling variance
df['Weights'] = 1 / df['Rolling_Var']

# Calculate weighted mean for detrended data
weighted_mean_detrended = np.average(df['Detrended'], weights=df['Weights'])

# Calculate weighted mean for original data
weighted_mean_original = np.average(df['Measurements'], weights=df['Weights'])

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(df['Time'], df['Measurements'], label='Original Measurements', alpha=0.6)
plt.plot(df['Time'], trend, label='Linear Trend', linestyle='--')
plt.plot(df['Time'], df['Detrended'], label='Detrended Measurements', alpha=0.6)
plt.xlabel('Time (minutes)')
plt.ylabel('Measurement Value')
plt.title('Time-Series Data with Linear Drift and Detrending')
plt.legend()
plt.show()

print(f"Weighted Mean of Detrended Data: {weighted_mean_detrended:.4f}")
print(f"Weighted Mean of Original Data: {weighted_mean_original:.4f}")
