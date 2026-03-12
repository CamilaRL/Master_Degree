import numpy as np
import matplotlib.pyplot as plt
from qutip import *


def Read_Density_Matrices(filename, dim):

    rhoTList = np.loadtxt(filename, dtype='complex')

    tlist = []
    rhot = []
    
    for t, rho in enumerate(rhoTList):
        
        tlist.append(t)
        rhot.append(Qobj(rho.reshape((dim, dim))))
    
    return tlist, rhot


def vonNeumann_Entropy(rho):
    
    S = []
    
    for rhot in rho:
        
        St = 0
        
        evals = Qobj(rhot).eigenenergies()
        
        for ev in evals:
            if ev > 0:
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


tempo_real = np.arange(0, 30, 0.01)

cnameList = ['min', 'max']

gList, cmin, cmax = np.loadtxt(f'./DensityMatrices/correlations.txt', unpack=True)

cList = [cmin, cmax]

for i, c in enumerate(cnameList):
    for j, g in enumerate(gList):

        ## reading rho(t)
        
        tempo_index, rho_q1_t = Read_Density_Matrices(f'./DensityMatrices/rhof_q1_c{c}_g{g}.txt', 2)
        tempo_index, rho_q2_t = Read_Density_Matrices(f'./DensityMatrices/rhof_q2_c{c}_g{g}.txt', 2)
        tempo_index, rho_total_t = Read_Density_Matrices(f'./DensityMatrices/rhof_c{c}_g{g}.txt', 4)
        
        
        MI = MutualInformation(rho_q1_t, rho_q2_t, rho_total_t)


        plt.plot(tempo_real, MI, label=f'|c| = {cList[i][j]:.3f} \n g = {g}')

    
plt.xlabel('Time')
plt.ylabel('Mutual Information')
plt.xscale('log')
plt.legend()
plt.show()


