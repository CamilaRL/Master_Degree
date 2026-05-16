import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os



def Trace_Distance(rho):

    D = []

    for t in range(len(rho)):
    
        drho = rho[t] - rho[-1]

        modrho = (drho.dag() * drho).sqrtm()
        
        D.append(modrho.tr())

    return D



## parameters

rho1_name = f'rhof_q1_t'
modo1 = 'Heating'

rho2_name = f'rhof_q2_t'
modo2 = 'Cooling'


tempo_real = np.arange(0, 20, 0.01)


## reading coherences
cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)


## reading time
tempo_index = []

for arquivo in os.listdir(f'./DensityMatrices/c_{cmod[0]}'):

    if rho1_name in arquivo:
    
        tempo_index.append(int(arquivo.replace(rho1_name, '').replace('.txt', '')))

tempo_index.sort()


for c in cmod:

    print(c)

    rho_q1_t_list = []
    rho_q2_t_list = []
    rho_total_t_list = []

    ## reading rho(t)
    
    for t in tempo_index:

        rho_q1_t = np.loadtxt(f'./DensityMatrices/c_{c}/{rho1_name}{t}.txt', dtype='complex', unpack=True)
        rho_q2_t = np.loadtxt(f'./DensityMatrices/c_{c}/{rho2_name}{t}.txt', dtype='complex', unpack=True)
        rho_total_t = np.loadtxt(f'./DensityMatrices/c_{c}/rhof_t{t}.txt', dtype='complex', unpack=True)
    
        rho_q1_t_list.append(Qobj(rho_q1_t))
        rho_q2_t_list.append(Qobj(rho_q2_t))
        rho_total_t_list.append(Qobj(rho_total_t))
        
    
    D_q1 = Trace_Distance(rho_q1_t_list)
    D_q2 = Trace_Distance(rho_q2_t_list)
    D_total = Trace_Distance(rho_total_t_list)


    plt.plot(tempo_real, D_total, color='black', label='Total')
    plt.plot(tempo_real, D_q1, color='red', label=modo1)
    plt.plot(tempo_real, D_q2, color='blue', label=modo2)
    
    plt.xlabel('Time')
    plt.ylabel('Trace Distance')
    plt.title(f'|c| = {c:.3f}')
    plt.xscale('log')
    plt.legend()
    plt.show()



