import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os


def vonNeumann_Entropy(rho):
    
    S = []
    
    for rhot in rho:
        
        St = 0
        
        evals = Qobj(rhot).eigenenergies()

        for ev in evals:
            St = St - ev * np.log(ev)
        
        S.append(St)

    return S
    

def MutualInformation(rho1, rho2, rho12):

    S1 = vonNeumann_Entropy(rho1)
    S2 = vonNeumann_Entropy(rho2)
    S12 = vonNeumann_Entropy(rho12)

    MI = []
    
    for t in range(len(S1)):
        
        MI.append(S1[t] + S2[t] - S12[t])

    return MI


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
        
    
    MI = MutualInformation(rho_q1_t_list, rho_q2_t_list, rho_total_t_list)


    plt.plot(tempo_real, MI, label=f'|c| = {c:.3f}')

    
plt.xlabel('Time')
plt.ylabel('Mutual Information')
plt.xscale('log')
plt.legend()
plt.show()


