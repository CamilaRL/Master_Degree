import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, kdq_r, kdq_i = np.loadtxt('./results/KDQ/kdq_zi_xi_00.txt', unpack=True)


## Create 2D scatter plot

plt.scatter(p, kdq_r, s=10)
plt.ylabel('KDQ')
plt.xlabel('p')
plt.title('Kirkwood-Dirac Quasiprobability for Werner States')
plt.show()

