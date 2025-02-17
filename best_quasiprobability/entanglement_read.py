import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

c1, c2, c3, S = np.loadtxt('./results/entanglement_entropy_10.txt', unpack=True)



fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

viridis_new = mpl.colormaps['viridis']
newcmap = ListedColormap(viridis_new(np.linspace(0, 0.1, 1)))

surf = ax.plot_surface(c1, c2, c3, cmap=newcmap, linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
