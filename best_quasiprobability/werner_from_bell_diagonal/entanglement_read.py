import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

c1, c2, c3, S = np.loadtxt('./results/entropy/entropy.txt', unpack=True)


## Create a 2D scatter plot

plt.scatter(c1, S, s=10)
plt.xlabel('a')
plt.ylabel(r'$E_f$')
plt.title('Formation Entropy')
plt.show()

## Create 3D scatter plot

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with color mapping
sc = ax.scatter(c1, c2, c3, c=S, cmap='plasma', marker='o')

# Labels and title
ax.set_xlabel("a1")
ax.set_ylabel("a2")
ax.set_zlabel("a3")
ax.set_title("Formation Entropy of Werner States")

# Add color bar
cbar = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
cbar.set_label("Formation Entropy")

plt.tight_layout()

# Show the plot
plt.show()

