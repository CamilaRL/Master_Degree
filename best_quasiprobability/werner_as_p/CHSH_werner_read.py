import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, chsh_r, chsh_i = np.loadtxt('./results/CHSH/CHSH_00.txt', unpack=True)


## Create a 2D scatter plot

plt.scatter(p, chsh_r, s=10, color='black')
plt.hlines(2, xmin=p[0], xmax=p[-1], color='red', label='CHSH limit')
plt.xlabel('p')
plt.ylabel(r'$CHSH$')
plt.title('CHSH for Werner States')
plt.legend()
plt.show()


