import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, qd_r, qd_i = np.loadtxt('./results/NDQP/ndqp_zz_xx_10.txt', unpack=True)

medidas = np.arange(0, 15, 1)

colors = plt.cm.viridis(np.linspace(0,1,len(p)))

## Create 2D scatter plot
for j in range(0, len(p)-15, 15):
	plt.plot(medidas, qd_r[j:j+15], color=colors[j], label=f'p = {p[j]:.2f}')

plt.xlabel('p')
plt.ylabel('NDQP')
plt.title('Nondemolition Quasiprobability for Werner States')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.show()

