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
ax.set_xlabel(r'$q$', fontsize=14)
ax.set_ylabel(r'$p$', fontsize=14)
ax.set_zlabel(r'$W(q,p)$', fontsize=14)


plt.show()
