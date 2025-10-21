import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os


def Quantum_Fisher_Information(rho_list, temp_list, dimensao):

    
    ### rho derivative respect to temperature
    
    size_rhoT = rho_list[0].shape[0]
    
    drho_list = [np.zeros((size_rhoT, size_rhoT), dtype=complex) for k in range(len(temp_list)-1)]
    
    for i in range(size_rhoT):
        for j in range(size_rhoT):
            
            for t in range(len(temp_list)-1):
    
                rho_i = rho_list[t].full()[i][j]
                rho_f = rho_list[t+1].full()[i][j]
                temp_i = temp_list[t]
                temp_f = temp_list[t+1]

                drho_list[t][i][j] = (rho_f - rho_i)/(temp_f - temp_i)
    
    ### Quantum Fisher Information calculation
    
    QFI_T = []
    
    for t in range(len(temp_list)-1):
        
        print(t)
        
        autoval, autovec = rho_list[t].eigenstates()
        
        QFI = complex(0,0)       
        
        for n in range(len(autoval)):
            for m in range(len(autoval)):
                
                bra = autovec[n].dag()
                ket = autovec[m]
                
                bra_drho_ket = (bra * Qobj(drho_list[t], dims=dimensao) * ket)

                QFI = QFI + (2*abs(bra_drho_ket)**2)/(abs(autoval[n]) + abs(autoval[m]))
        
        QFI_T.append(QFI.real)
        
    return temp_list[:-1], QFI_T
    

  
## read rhot files

tempo = []

for arquivo in os.listdir('./DensityMatrices/'):

    tempo.append(int(arquivo.replace('rhof_t', '').replace('.txt', '')))

tempo.sort()

rhot_list = []

for t in tempo:

    rhot = np.loadtxt(f'./DensityMatrices/rhof_t{t}.txt', usecols=(0,1,2,3), dtype='complex')

    rhot_list.append(Qobj(rhot))


## Compute QFI

tempo_list = np.arange(0.005, 5, 0.001)

tempo, QFI = Quantum_Fisher_Information(rhot_list, tempo_list, [[4],[4]])


## save QFI in file

f = open('./FisherInformation/QFI.txt', 'w')

for t in range(len(tempo)):

    f.write(f'{tempo[t]} {QFI[t]}\n')
    
    
f.close()










