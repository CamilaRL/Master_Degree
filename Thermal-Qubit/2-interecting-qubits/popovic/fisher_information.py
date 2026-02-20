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
        
        autoval, autovec = rho_list[t].eigenstates()
        
        QFI = 0 + 0j      
        
        for n in range(len(autoval)):
            for m in range(len(autoval)):
                
                bra = autovec[n].dag()
                ket = autovec[m]
                
                bra_drho_ket = bra * Qobj(drho_list[t], dims=dimensao) * ket
                
                QFI = QFI + (2*abs(bra_drho_ket)**2)/(abs(autoval[n]) + abs(autoval[m]))
                
        QFI_T.append(QFI.real)
        
    return temp_list[:-1], QFI_T
    

  
## parameters

qubit = 'q1'
modo = 'Heating'

os.mkdir(f'./FisherInformation/{modo}')

tempo_real = np.arange(0, 30, 0.01)


## reading density matrices

JxList = np.loadtxt('./DensityMatrices/T_final_qubit.txt', unpack=True, usecols=(0))

for Jx in JxList:

    tempo_index, rhof_q_min = Read_Density_Matrices(f'./DensityMatrices/rhof_{qubit}_cmin_Jx{Jx}.txt', 2)
    tempo_index, rhof_q_max = Read_Density_Matrices(f'./DensityMatrices/rhof_{qubit}_cmax_Jx{Jx}.txt', 2)

    ## Compute QFI

    tempo_QFI_min, QFI_min = Quantum_Fisher_Information(rhof_q_min, tempo_real, [[2],[2]])
    tempo_QFI_max, QFI_max = Quantum_Fisher_Information(rhof_q_max, tempo_real, [[2],[2]])

    ## save QFI in file

    fmin = open(f'./FisherInformation/{modo}/QFI_{qubit}_cmin_Jx{Jx}.txt', 'w')
    fmax = open(f'./FisherInformation/{modo}/QFI_{qubit}_cmax_Jx{Jx}.txt', 'w')

    for t in range(len(tempo_QFI_min)):

        fmin.write(f'{tempo_QFI_min[t]} {QFI_min[t]}\n')
        fmax.write(f'{tempo_QFI_max[t]} {QFI_max[t]}\n')
        
    fmin.close()
    fmax.close()










