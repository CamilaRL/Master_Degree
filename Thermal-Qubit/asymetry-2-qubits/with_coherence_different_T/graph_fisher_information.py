import numpy as np
import matplotlib.pyplot as plt

qubit = '2'
modo = 'Cooling'

cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, len(cmod))))

cmod.sort()

for c in cmod:

    tempo, QFI = np.loadtxt(f'./FisherInformation/{modo}/QFI_q{qubit}_c{c}.txt', unpack=True)

    cor = next(colors)
    
    plt.plot(tempo, QFI, color=cor,  label=f'{c}')


plt.xscale('log')
plt.legend(loc='upper right')
plt.ylabel('Quantum Fisher Information')
plt.xlabel('Time')
plt.title(f'{modo} - Qubit {qubit}')
plt.tight_layout()
plt.show()
