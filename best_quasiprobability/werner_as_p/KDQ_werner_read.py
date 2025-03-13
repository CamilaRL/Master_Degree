import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, kdq_r, kdq_i = np.loadtxt('./results/KDQ/kdq_zz_xx_00.txt', unpack=True)

medidas = np.arange(0, 15, 1)

colors = plt.cm.viridis(np.linspace(0,1,len(p)))

for j in range(0, len(p)-15, 15):

	plt.plot(medidas, kdq_r[j:j+15], color=colors[j], label=f'p = {p[j]:.2f}')
	

plt.xlabel("Projetores")
plt.ylabel("KDQ Real")

plt.suptitle("KDQ Distribution of Werner States")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()
plt.show()

'''
# Create 3D scatter plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

for j in range(0, len(p)-15, 15):
# Scatter plot with color mapping
	sc = ax.scatter(projetores, kdq_r[j:j+15], kdq_i[j:j+15], marker='o', label=f'p = {p[j]:.2f}')

# Labels and title
ax.set_xlabel("Projetores")
ax.set_ylabel("KDQ Real")
ax.set_zlabel("KDQ Imaginary")
ax.set_title("KDQ Distribution of Werner States")

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()

# Show the plot
plt.show()
'''
