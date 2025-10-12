import numpy as np
import matplotlib.pyplot as plt


modoList = ['Cooling', 'Heating']
dSr = 0.1

colors = ['blue', 'red']

for m, modo in enumerate(modoList):
    curvas = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
    cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)
    
    for i in range(len(curvas)):
        
        tlist, Srt = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/Sr_curve_{i}.txt', unpack=True)
        
        plt.plot(tlist, Srt, color=colors[m], label=f'{modo} - |c| = {cmod[i]:.1f}')


plt.legend(loc='center right')
plt.ylabel('Relative Entropy')
plt.xlabel('Time')
plt.title(modo)
plt.show()
    
    
