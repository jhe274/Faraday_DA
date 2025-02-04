import numpy as np
import matplotlib.pyplot as plt

# Generate sample data (asymmetric envelope distribution)
np.random.seed(42)
data = np.concatenate([
    np.random.normal(loc=50, scale=10, size=700),  # Main peak
    np.random.normal(loc=30, scale=8, size=300)    # Smaller secondary peak (asymmetry)
])

# Define histogram bins
num_bins = 20
counts, bin_edges = np.histogram(data, bins=num_bins)

# Define the baseline function (linear)
slope = 0.1  # Adjust as needed
offset = 5    # Adjust as needed
baseline = slope * bin_edges[:-1] + offset  # Compute baseline values

# Compute asymmetric upper and lower extensions
median_value = np.median(counts)  # Use median to split upper/lower variations
upper_extension = np.maximum(counts - median_value, 0)  # How much the data extends above the median
lower_extension = np.maximum(median_value - counts, 0)  # How much the data extends below the median

# Plot the histogram with an asymmetric envelope
fig, ax = plt.subplots(figsize=(8, 5))

# Draw histogram bars extending asymmetrically from the baseline
for i in range(len(counts)):
    x_left = bin_edges[i]
    x_right = bin_edges[i+1]

    # Upper and lower bounds for asymmetric envelope
    y_lower = baseline[i] - lower_extension[i]
    y_upper = baseline[i] + upper_extension[i]

    # Fill the asymmetric histogram
    ax.fill_between([x_left, x_right], y_lower, y_upper, color='b', alpha=0.6)

# Plot the baseline
ax.plot(bin_edges[:-1], baseline, 'r--', label="Baseline (y = {:.2f}x + {:.2f})".format(slope, offset))

# Labels and title
ax.set_xlabel("Value")
ax.set_ylabel("Frequency")
ax.set_title("Asymmetric Histogram with Envelope Baseline")
ax.legend()
plt.show()