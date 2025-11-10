import numpy as np
import matplotlib.pyplot as plt

qubit = '1'
modo = 'Heating'

tempo, QFI = np.loadtxt(f'./FisherInformation_{modo}/QFI_q{qubit}_c0.txt', unpack=True)

plt.plot(tempo, QFI)
plt.ylabel('Quantum Fisher Information')
plt.xlabel('Time')
plt.title(f'{modo} - Qubit {qubit}')
plt.show()
