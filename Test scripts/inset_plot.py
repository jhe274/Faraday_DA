import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

# Generate synthetic sine wave data
x = np.linspace(0, 2*np.pi, 1000)  # 0 to 2π with 1000 points
y = np.sin(x)

# Create the main figure and axes
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the sine wave
ax.plot(x, y, label=r"$y = \sin(x)$", color='b')

# Labels and title
ax.set_xlabel("x (radians)", fontsize=14)
ax.set_ylabel("y", fontsize=14)
ax.set_title("Sine Wave with Zoomed-In First Quarter", fontsize=16)
ax.legend()

# Create an inset plot at the lower left and zoom into the first quarter-period
axins = zoomed_inset_axes(ax, zoom=3, loc="lower left")  # 3x zoom
axins.plot(x, y, color='b')

# Define the zoom-in region (first quarter-period)
x1, x2 = 0, np.pi / 4  # First quarter: [0, π/2]
y1, y2 = 0, 0.1  # Zoom in on positive y-values

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Set ticks and hide labels for the inset
axins.xaxis.get_major_locator().set_params(nbins=5)
axins.yaxis.get_major_locator().set_params(nbins=5)
axins.tick_params(labelleft=False, labelbottom=False)

# Mark the zoomed region with a dashed rectangle
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="black", linestyle="dashed")

plt.show()
