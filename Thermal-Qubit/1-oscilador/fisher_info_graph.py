import numpy as np
import matplotlib.pyplot as plt


tot = 13
cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, tot)))

curvas = []
completion1 = []
completion2 = []

for i in range(1, tot, 1):

    tlist, QFI = np.loadtxt(f'./FisherInformation_Aquecer/curva_{i}.txt', unpack=True)
    
    cor = next(colors)
    
    plt.plot(tlist, QFI, color=cor, label=f'Curva {i}')
    plt.xlabel('Time')
    plt.ylabel('Quantum Fisher Information')

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()
