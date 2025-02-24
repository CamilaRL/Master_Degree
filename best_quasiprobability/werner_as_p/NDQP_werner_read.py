import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, qd_r, qd_i = np.loadtxt('./results/NDQP/ndqp_zz_xx_00.txt', unpack=True)


## Create 2D scatter plot

plt.scatter(p, qd_r, s=10)
plt.xlabel('p')
plt.ylabel('NDQP')
plt.title('Nondemolition Quasiprobability for Werner States')
plt.show()

