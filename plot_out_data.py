import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('cmake-build-debug/out_data.csv', delimiter=',', skiprows=1)

plt.plot(data[:, 0], data[:, 1])
plt.show()
