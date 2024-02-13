import numpy as np
import os
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt 


cwd = os.getcwd()
data = np.loadtxt(f'{os.getcwd()}/data.txt')
# print(data)
data = data[np.argsort(data[:, 0])]
print(data)

x, y = data.T
ys_cs = CubicSpline(x, y)
xx = np.linspace(x.min(), x.max(), 1000)
yy = ys_cs(xx)

savgol_params = {'window_length': 11, 'polyorder': 3}

ys = savgol_filter(y, **savgol_params)
dy = np.gradient(y, x)
dyy = np.gradient(yy, xx)
dyys = savgol_filter(yy, deriv=1, delta=xx[1]-xx[0], **savgol_params)
print(np.diff(x))
print(y.shape)

fig, axs = plt.subplots(2, 1, figsize=(12, 12), dpi=100)
axs[0].plot(x, y)
axs[0].plot(x, ys, '-.')
axs[0].plot(xx, yy, '-o')

# axs[1].plot(x, dy)
axs[1].plot(xx, dyy, '-o')
axs[1].plot(xx, dyys)
axs[0].tick_params(axis='x', labelsize=15)
axs[0].tick_params(axis='y', labelsize=15)
axs[0].set_ylabel('Intensity (au)', fontsize=20)
axs[0].set_xlabel('Frequency (GHz)', fontsize=20)
plt.show()