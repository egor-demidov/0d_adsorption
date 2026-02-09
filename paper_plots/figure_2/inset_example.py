import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Sample data
x = np.linspace(0, 10, 500)
y = np.sin(x) * np.exp(-0.1 * x)

fig, ax = plt.subplots()

# Main plot
ax.plot(x, y)
ax.set_title("Main Plot with Zoomed Inset")

# Create inset axes
axins = inset_axes(ax, width="40%", height="30%", loc="upper right")

# Plot same data in inset
axins.plot(x, y)

# Define zoom region
x1, x2 = 2, 4
y1, y2 = -0.5, 0.5
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Remove ticks (optional)
axins.set_xticks([])
axins.set_yticks([])

# Draw box and connectors
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.show()
