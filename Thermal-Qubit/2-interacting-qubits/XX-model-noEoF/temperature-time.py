import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os


def Read_Density_Matrices(filename, dim):

    rhoTList = np.loadtxt(filename, dtype='complex')

    tlist = []
    rhot = []
    
    for t, rho in enumerate(rhoTList):
        
        tlist.append(t)
        rhot.append(Qobj(rho.reshape((dim, dim))))
    
    return tlist, rhot




    
## parameters

tempo_real = np.arange(0, 30, 0.01)

w0 = 1

## reading coherences
gList = [0.0, 0.8]
cList = ['cmin', 'cmax']
modoList = ['Heating', 'Cooling']

temp_q1 = []
temp_q2 = []

for c in cList:

    for g in gList:
    ## reading rho(t)
    
        temp_aux1 = []
        temp_aux2 = []
    
        tempo_index, rho_q1_t = Read_Density_Matrices(f'./DensityMatrices/rhof_q1_{c}_g{g}.txt', 2)
        tempo_index, rho_q2_t = Read_Density_Matrices(f'./DensityMatrices/rhof_q2_{c}_g{g}.txt', 2)
        
        for t in tempo_index:
            
            beta_q1 = (1/w0)*np.log(rho_q1_t[t].full()[0][0].real/rho_q1_t[t].full()[1][1].real)
            
            beta_q2 = (1/w0)*np.log(rho_q2_t[t].full()[0][0].real/rho_q2_t[t].full()[1][1].real)
            
            temp_aux1.append(beta_q1)
            temp_aux2.append(beta_q2)
        
        temp_q1.append(temp_aux1)
        temp_q2.append(temp_aux2)
    

fig = plt.figure(figsize=(10,5))

plt.subplot(121)
    
plt.plot(tempo_real, temp_q1[0], color='red', linestyle='--', linewidth=2, label='J = 0 - Heating')
plt.plot(tempo_real, temp_q2[0], color='blue', linestyle='--', linewidth=2, label='J = 0 - Cooling')
    
plt.plot(tempo_real, temp_q1[1], color='red', linewidth=2, label='J = 0.8 - Heating')
plt.plot(tempo_real, temp_q2[1], color='blue', linewidth=2, label='J = 0.8 - Cooling')

plt.xlabel('Time', fontsize=14)
plt.ylabel('Qubit\'s Temperature', fontsize=14)
plt.title('Zero Initial Coherence', fontsize=16)
plt.legend(loc='upper right', fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)


plt.subplot(122)
    
plt.plot(tempo_real, temp_q1[2], color='red', linestyle='--', linewidth=2, label='J = 0 - Heating')
plt.plot(tempo_real, temp_q2[2], color='blue', linestyle='--', linewidth=2, label='J = 0 - Cooling')
    
plt.plot(tempo_real, temp_q1[3], color='red', linewidth=2, label='J = 0.8 - Heating')
plt.plot(tempo_real, temp_q2[3], color='blue', linewidth=2, label='J = 0.8 - Cooling')

plt.xlabel('Time', fontsize=14)
plt.title('Maximum Initial Coherence', fontsize=16)
plt.legend(loc='upper right', fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)


plt.tight_layout()
plt.show()