import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os
import math


def Read_Density_Matrices(filename, dim):

    rhoTList = np.loadtxt(filename, dtype='complex')

    tlist = []
    rhot = []
    
    for t, rho in enumerate(rhoTList):
        
        tlist.append(t)
        rhot.append(Qobj(rho.reshape((dim, dim)), dims=[[2,2],[2,2]]))
    
    return tlist, rhot


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

J = 1
delta = -2

tempo_real = np.arange(0, 30, 0.01)


cList = ['cmin', 'cmax']
cName = ['Zero Initial Coherence', 'Maximum Initial Coherence']

gList = np.loadtxt(f'./DensityMatrices_J{J}_delta{delta}/T_final_qubit.txt', unpack=True, usecols=(0))

cmap = plt.get_cmap('rainbow')

for c in range(2):

    colors = iter(cmap(np.linspace(0.1, 1, len(gList))))
    
    for g in gList:

        tempo_index, rho_total = Read_Density_Matrices(f'./DensityMatrices_J{J}_delta{delta}/rhof_{cList[c]}_g{g}.txt', 4)
    
        E = EoF(rho_total[:-1])

        cor = next(colors)
        plt.plot(tempo_real, E, color=cor, label=f'g = {g}')

    plt.xlabel('Time')
    plt.ylabel('Entanglement of Formation')
    plt.title(cName[c])
    plt.legend()
    plt.show()


