import numpy as np
import matplotlib.pyplot as plt
import os



### MAIN ###

modo = 'Heating'
dSr = 0.1

curvas, cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True)

t_equilibrio = []

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(curvas))))

for i, j in enumerate(curvas):

    curva = int(j)
    
    tlist, completion = np.loadtxt(f'./ThermalKinematics_{modo}_{dSr}/completion_{curva}.txt', unpack=True)
    
    equilibrio = np.where(completion==completion[-1])
    
    t_equilibrio.append(tlist[equilibrio[0][0]])
    
    cor = next(colors)
    
    plt.plot(tlist, completion, color=cor, label=f'|c| = {cmod[i]:.3f}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.title(modo)
plt.xscale('log')
plt.xlim(left=0.01)
plt.tight_layout()
plt.show()

