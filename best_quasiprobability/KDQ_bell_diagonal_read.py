import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

a1, a2, a3, kdq_r, kdq_i = np.loadtxt('./results/KDQ/kdq_zi_yi.txt', unpack=True)


# Create 3D scatter plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with color mapping
sc = ax.scatter(a1, a2, a3, c=kdq_r, cmap='plasma', marker='o')

# Labels and title
ax.set_xlabel("a1")
ax.set_ylabel("a2")
ax.set_zlabel("a3")
ax.set_title("KDQ Distribution of Bell Diagonal States")

# Add color bar
cbar = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
cbar.set_label("KDQ")

plt.tight_layout()

# Show the plot
plt.show()

