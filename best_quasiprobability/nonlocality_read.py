import matplotlib.pyplot as plt
import numpy as np

a1, a2, HSD = np.loadtxt('./results/nonlocality.txt', unpack=True)

A1, A2 = np.meshgrid(a1, a2)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(A1, A2, HSD)

plt.show()
