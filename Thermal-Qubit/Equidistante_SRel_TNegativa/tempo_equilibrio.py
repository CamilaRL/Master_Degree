import numpy as np
import matplotlib.pyplot as plt
import os

    

### MAIN ###

Sr = 0.1
modoList = ['Cooling', 'Heating']
cores = ['blue', 'red']


for modo in modoList:
    
    m = modoList.index(modo)
    
    curvas = np.loadtxt(f'./p_final_0.2/FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
    cmod = np.loadtxt(f'./p_final_0.2/FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)
    
    t_equilibrio = []
    
    cmap = plt.get_cmap('rainbow')
    colors = iter(cmap(np.linspace(0.01, 1, len(curvas))))

    for i in range(len(curvas)):

        tlist, completion = np.loadtxt(f'./p_final_0.2/ThermalKinematics_{modo}_{Sr}/completion_{int(curvas[i])}.txt', unpack=True)        
        
        equilibrio = np.where(completion==completion[-1])
        
        t_equilibrio.append(tlist[equilibrio[0][0]])
        

    plt.scatter(cmod, t_equilibrio, color=cores[m])
    plt.plot(cmod, t_equilibrio, color=cores[m], label=f'{modoList[m]}')

plt.ylabel('Equilibration Time')
plt.xlabel('|c|')
plt.legend(loc='center right')
plt.tight_layout()
plt.show()
