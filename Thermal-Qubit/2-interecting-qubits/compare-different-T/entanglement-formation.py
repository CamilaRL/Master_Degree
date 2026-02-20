import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os
import math



def EoF(rho):

    E = []

    for t in range(len(rho)):
        
        C = concurrence(rho[t])
        
        x = 0.5*(1 + np.sqrt(1 - C**2))

        if x == 0 or x == 1:
            h = 0
            
        else:
            h = -x * math.log2(x) - (1-x) * math.log2(1-x)
        
        E.append(h)
        
    return E
        
    

## parameters

rho1_name = f'rhof_q1_t'
modo1 = 'Heating'

rho2_name = f'rhof_q2_t'
modo2 = 'Cooling'


tempo_real = np.arange(0, 20, 0.01)


## reading coherences
cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)
cmod.sort()


## reading time
tempo_index = []

for arquivo in os.listdir(f'./DensityMatrices/c_{cmod[0]}'):

    if rho1_name in arquivo:
    
        tempo_index.append(int(arquivo.replace(rho1_name, '').replace('.txt', '')))

tempo_index.sort()


for c in cmod:

    print(c)

    rho_total_t_list = []

    ## reading rho(t)
    
    for t in tempo_index:

        rho_total_t = np.loadtxt(f'./DensityMatrices/c_{c}/rhof_t{t}.txt', dtype='complex', unpack=True)
    
        rho_total_t_list.append(Qobj(rho_total_t, dims=[[2, 2], [2, 2]]))
        
    
    E = EoF(rho_total_t_list)


    plt.plot(tempo_real, E, label=f'|c| = {c:.3f}')

    
plt.xlabel('Time')
plt.ylabel('Entanglement of Formation')
plt.xscale('log')
plt.legend()
plt.show()

