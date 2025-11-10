import numpy as np
import matplotlib.pyplot as plt
import os
    

### MAIN ###

dSr = 0.1
modo = 'Heating'

curvas = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(curvas))))


for i, j in enumerate(curvas):

    curva = int(j)    
    
    tlist, QFI = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/QFI_curve_{curva}.txt', unpack=True)
    
    cor = next(colors)
    
    plt.plot(tlist, QFI, color=cor, label=f'|c| = {cmod[i]:.3f}')
    plt.xlabel('Time')
    plt.ylabel('Quantum Fisher Information')

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.title('Heating')
plt.tight_layout()
plt.show()
