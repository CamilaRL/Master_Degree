import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, S = np.loadtxt('./results/entropy/entropy_00.txt', unpack=True)


## Create a 2D scatter plot

plt.scatter(p, S, s=10)
plt.xlabel('p')
plt.ylabel(r'$E_f$')
plt.title('Formation Entropy for Werner States')
plt.show()


