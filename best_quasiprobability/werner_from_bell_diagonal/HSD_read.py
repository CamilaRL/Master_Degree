import matplotlib.pyplot as plt
import numpy as np

a1, a2, a3, HSD = np.loadtxt('./results/HSD/HSD.txt', unpack=True)


# Create 2D scatter plot

plt.scatter(a1, HSD, s=10)
plt.ylabel('HSD')
plt.xlabel('a')
plt.title('Hilber-Schmidt Distance for Werner States')
plt.show()


# Create 3D scatter plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with color mapping
sc = ax.scatter(a1, a2, a3, c=HSD, cmap='plasma', marker='o')

# Labels and title
ax.set_xlabel("a1")
ax.set_ylabel("a2")
ax.set_zlabel("a3")
ax.set_title("Hilbert-Schmidt Distace of Bell Diagonal States")

# Add color bar
cbar = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
cbar.set_label("HSD")

plt.tight_layout()

# Show the plot
plt.show()

