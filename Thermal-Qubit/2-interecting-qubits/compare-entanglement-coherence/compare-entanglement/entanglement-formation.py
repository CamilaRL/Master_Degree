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

    for t in [0]:#range(len(rho)):
        rhot = rho[t]
        
        lambda1 = abs(rhot[2][1]) - np.sqrt(rhot[0][0]*rhot[3][3])
        lambda2 = abs(rhot[3][0]) - rhot[1][1]
        print(lambda1, lambda2)
        '''rhot = Qobj(rho[t], dims=[[2, 2], [2, 2]])
        
        sysy = tensor(sigmay(), sigmay())

        rho_tilde = (rhot * sysy) * (rhot.conj() * sysy)

        evals = rho_tilde.eigenenergies()
        print(evals)
        # abs to avoid problems with sqrt for very small negative numbers
        evals = abs(np.sort(np.real(evals)))
        print(evals)
        lsum = np.sqrt(evals[3]) - np.sqrt(evals[2]) - np.sqrt(evals[1]) - np.sqrt(evals[0])

        print(max(0, lsum))
        
        #C = concurrence(rhot)
        #print(C)
        x = 0.5*(1 + np.sqrt(1 - C**2))

        if x == 0 or x == 1:
            h = 0
            
        else:
            h = -x * math.log2(x) - (1-x) * math.log2(1-x)
        
        E.append(h)
        
    return E'''
        
    

## parameters

tempo_real = np.arange(0, 30, 0.01)

cnameList = ['min', 'max']

gList, cmin, cmax = np.loadtxt(f'./DensityMatrices/correlations.txt', unpack=True)

cList = [cmin, cmax]

for i, c in enumerate(cnameList):
    for j, g in enumerate(gList):
    
        ## reading rho(t)
        tempo_index, rho_total_t = Read_Density_Matrices(f'./DensityMatrices/rhof_c{c}_g{g}.txt', 4)

        E = EoF(rho_total_t)

        '''plt.plot(tempo_real, E, label=f'|c| = {cList[i][j]:.3f} \n g = {g}')

    
plt.xlabel('Time')
plt.ylabel('Entanglement of Formation')
plt.xscale('log')
plt.legend()
plt.show()
'''
