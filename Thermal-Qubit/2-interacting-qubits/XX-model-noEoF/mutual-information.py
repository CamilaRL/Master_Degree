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

colors_coherence = ['teal', 'purple']

lines_coupling = ['--', '-']

gList, cmin, cmax = np.loadtxt(f'./DensityMatrices/correlations.txt', unpack=True)

titulo = ['Zero Initial Coherence', 'Maximum Initial Coherence']

cList = [cmin, cmax]

handles = []
fig = plt.figure()
for i, c in enumerate(cnameList):
    color_handle = []
    for j, g in enumerate(gList):

        ## reading rho(t)
        
        tempo_index, rho_q1_t = Read_Density_Matrices(f'./DensityMatrices/rhof_q1_c{c}_g{g}.txt', 2)
        tempo_index, rho_q2_t = Read_Density_Matrices(f'./DensityMatrices/rhof_q2_c{c}_g{g}.txt', 2)
        tempo_index, rho_total_t = Read_Density_Matrices(f'./DensityMatrices/rhof_c{c}_g{g}.txt', 4)
        
        
        MI = MutualInformation(rho_q1_t, rho_q2_t, rho_total_t)


        mi_plot = plt.plot(tempo_real[:2000], MI[:2000], color=colors_coherence[i], linestyle=lines_coupling[j], linewidth=2, label=f'{titulo[i]} \n J = {g}')
        color_handle.append(mi_plot[0])
        
    handles.append(color_handle) 

labels = [f'J = {gList[0]:.1f}', f'J = {gList[1]:.1f}']

# Legenda Heating (Vermelha) - Superior
leg_h = fig.legend(handles[0], labels, 
                   loc='lower center', 
                   ncol=3, 
                   title=f'{titulo[0]}', 
                   title_fontproperties={'weight':'bold', 'size':12},
                   fontsize=12,
                   bbox_to_anchor=(0.5, 0.1), 
                   frameon=False)

# Legenda Cooling (Azul) - Inferior
leg_c = fig.legend(handles[1], labels, 
                   loc='lower center', 
                   ncol=3, 
                   title=f'{titulo[1]}', 
                   title_fontproperties={'weight':'bold', 'size':12},
                   fontsize=12,
                   bbox_to_anchor=(0.5, 0.01), 
                   frameon=False)
    
plt.xlabel('Time', fontsize=12)
plt.ylabel('Mutual Information', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.tight_layout()
plt.subplots_adjust(bottom=0.35)
plt.show()


