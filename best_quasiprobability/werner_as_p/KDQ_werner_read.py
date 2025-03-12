import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, kdq_r, kdq_i = np.loadtxt('./results/KDQ/kdq_zz_xx_11.txt', unpack=True)


# Create 3D scatter plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with color mapping
sc = ax.scatter(p, kdq_r, kdq_i, marker='o')

# Labels and title
ax.set_xlabel("p")
ax.set_ylabel("KDQ Real")
ax.set_zlabel("KDQ Imaginary")
ax.set_title("KDQ Distribution of Werner States")

plt.tight_layout()

# Show the plot
plt.show()
