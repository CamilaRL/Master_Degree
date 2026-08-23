import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

x=np.linspace(-10,10, num=100)
y=np.linspace(-10,10, num=100)

x, y = np.meshgrid(x, y)

z = np.exp(-0.1*x**2-0.1*y**2)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x,y,z, cmap=cm.jet)
ax.set_xlabel(r'$\alpha$', fontsize=14)
ax.set_ylabel(r'$\alpha^*$', fontsize=14)

ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$W(\mathbf{x})$', fontsize=14, rotation=0, labelpad=10)


plt.show()
