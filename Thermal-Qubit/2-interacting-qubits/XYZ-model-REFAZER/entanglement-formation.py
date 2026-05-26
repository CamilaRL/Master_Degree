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
        rhot.append(Qobj(rho.reshape((dim, dim))))
    
    return tlist, rhot


def EoF(rho):

    E = []

    for t in range(len(rho)):
        rhot = Qobj(rho[t], dims=[[2, 2], [2, 2]])
        
        C = concurrence(rhot)
        x = 0.5*(1 + np.sqrt(1 - C**2))

        if x == 0 or x == 1:
            h = 0
            
        else:
            h = -x * math.log2(x) - (1-x) * math.log2(1-x)
        
        E.append(h)
        
    return E
        
    

## parameters

tempo_real = np.arange(0, 30, 0.01)

cnameList = ['min', 'max']

colors_coherence = ['teal', 'purple']

lines_coupling = ['-', '--']

gList, cmin, cmax = np.loadtxt(f'./DensityMatrices/correlations.txt', unpack=True)

cList = [cmin, cmax]

for j, g in enumerate(gList):

    for i, c in enumerate(cnameList):
    
        ## reading rho(t)
        tempo_index, rho_total_t = Read_Density_Matrices(f'./DensityMatrices/rhof_c{c}_g{g}.txt', 4)
        
        #rhof = rho_total_t[-1].full()
        
        #f rhof[1][1]*rhof[2][2] < abs(rhof[0][3])**2 or rhof[0][0]*rhof[3][3] < abs(rhof[1][2])**2:
        #    print(':) Emaranhado!')
        #else:
        #    print(':( Não emaranhado')
        
        E = EoF(rho_total_t)

        plt.plot(tempo_real, E, color=colors_coherence[i], linestyle=lines_coupling[i], linewidth=2, label=f'|c| = {cList[i][j]:.3f}')
    
    
    plt.title(f'J = {g}')
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Entanglement of Formation', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()
