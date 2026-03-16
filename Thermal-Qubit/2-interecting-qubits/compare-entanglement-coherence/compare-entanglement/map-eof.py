import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os
import math


def EoF(rho):

    rho = Qobj(rho, dims=[[2, 2], [2, 2]])
        
    C = concurrence(rho)

    x = 0.5*(1 + np.sqrt(1 - C**2))

    if x == 0 or x == 1:
        E = 0
            
    else:
        E = -x * math.log2(x) - (1-x) * math.log2(1-x)
        
    return E
    

def Rho(w0, beta_1, beta_2, e, f):
    
    Z1 = 2*np.cosh(beta_1*w0/2)
    Z2 = 2*np.cosh(beta_2*w0/2)

    p1 = np.exp(beta_1*w0/2)/Z1
    p2 = np.exp(beta_2*w0/2)/Z2

    rho0_q1 = Qobj([[p1, 0],[0, 1-p1]])

    rho0_q2 = Qobj([[p2, 0],[0, 1-p2]])

    G = basis(2,0)
    E = basis(2,1)

    v00 = tensor(G, G)
    v01 = tensor(G, E)
    v10 = tensor(E, G)
    v11 = tensor(E, E)

    coherences_matrix = e * v01 * v10.dag() + e.conjugate() * v10 * v01.dag()
    entanglement_matrix = f * v00 * v11.dag() + f.conjugate() * v11 * v00.dag()

    rho0 = tensor(rho0_q1, rho0_q2) + coherences_matrix + entanglement_matrix
    
    return rho0
    



print(tensor(sigmay(), sigmay()))

### MAIN ###
'''
gList, cmin, cmax = np.loadtxt('./DensityMatrices/correlations.txt', unpack=True)
gList, Tqubit_final, Th, Tc = np.loadtxt('./DensityMatrices/temperature_qubit.txt', unpack=True)


for k, g in enumerate(gList):
    
    cList = np.arange(cmin[k], cmax[k], 0.001)
    cLabel = [cList[i] for i in range(0, len(cList), 15)]
    
    EoF_map = np.zeros((len(cList), len(cList)))
    
    for i, e in enumerate(cList):
        for j, f in enumerate(cList):
            
            rho0 = Rho(1, 1/Tc[k], 1/Th[k], e, f)
            if EoF(rho0) != 0:
                
                print('Emaranhamento!!!')
            EoF_map[i][j] = EoF(rho0)
            
    
    
    image = plt.imshow(EoF_map, cmap='viridis', aspect='equal', origin='lower', vmin=0)

    plt.yticks(ticks=np.arange(0, len(cList), 15), labels=cLabel)
    plt.xticks(ticks=np.arange(0, len(cList), 15), labels=cLabel)
    plt.colorbar()
    plt.title(f'Entanglement of Formation | g = {g}')
    plt.ylabel('e')
    plt.xlabel('f')
    plt.tight_layout()
    plt.show()'''