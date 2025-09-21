import numpy as np
import matplotlib.pyplot as plt


modo = 'Heating'
dSr = 0.1

curvas, cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmod))))    


for i in range(len(curvas)):
    
    tlist, Srt = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/Sr_curve_{i}.txt', unpack=True)
    
    cor = next(colors)
    
    plt.plot(tlist, Srt, color=cor, label=f'|c| = {cmod[i]:.4f}')


plt.legend(loc='center right')
plt.ylabel('Relative Entropy')
plt.xlabel('Time')
plt.title(modo)
plt.show()
    
    