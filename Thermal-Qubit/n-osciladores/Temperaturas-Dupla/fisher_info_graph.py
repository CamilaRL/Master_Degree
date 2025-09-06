import numpy as np
import matplotlib.pyplot as plt
import os
    

### MAIN ###

modo = 'Aquecer'

curvas, cmod = np.loadtxt(f'./FisherInformation_{modo}/cmod.txt', unpack=True)
cmod, curvas = (list(t) for t in zip(*sorted(zip(cmod, curvas))))

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(curvas))))


for i, j in enumerate(curvas):

    curva = int(j)    
    
    tlist, QFI = np.loadtxt(f'./FisherInformation_{modo}/QFI_curve_{curva}.txt', unpack=True)
    
    cor = next(colors)
    
    plt.plot(tlist, QFI, color=cor, label=f'|c| = {cmod[i]:.3f}')
    plt.xlabel('Time')
    plt.ylabel('Quantum Fisher Information')

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.title('Heating')
plt.tight_layout()
plt.show()
