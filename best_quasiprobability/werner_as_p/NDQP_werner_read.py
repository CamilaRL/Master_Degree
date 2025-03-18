import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

p, qd_r, qd_i = np.loadtxt('./results/NDQP/ndqp_zz_xx_00.txt', unpack=True)

medidas = np.arange(0, 15, 1)

colors = plt.cm.viridis(np.linspace(0,1,len(p)))

qd_r_mean = []
qd_i_mean = []
p_list = []

## Create 2D scatter plot
for j in range(0, len(p)-15, 15):
	#plt.plot(medidas, qd_r[j:j+15], color=colors[j], label=f'p = {p[j]:.2f}')
	
	qd_r_mean.append(sum(qd_r[j:j+15])/15)
	qd_i_mean.append(sum(qd_i[j:j+15])/15)
	p_list.append(p[j+15])

'''
plt.xlabel('p')
plt.ylabel('NDQP')
plt.title('Nondemolition Quasiprobability for Werner States')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.show()
'''

plt.scatter(p_list, qd_r_mean, s=10, label='Real')
plt.plot(p_list, qd_r_mean)
plt.scatter(p_list, qd_i_mean, s=10, label='Imaginary')
plt.plot(p_list, qd_i_mean)
plt.title("NDQP of Werner States")
plt.xlabel('p')
plt.ylabel('Mean NDQP')
plt.legend()
plt.show()
