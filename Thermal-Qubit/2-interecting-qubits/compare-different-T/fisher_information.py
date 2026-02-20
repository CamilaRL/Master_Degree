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
rho_name = f'rhof_{qubit}_t'
modo = 'Heating'

os.mkdir(f'./FisherInformation/{modo}')

tempo_real = np.arange(0, 20, 0.01)


## reading coherences
cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)


## reading time
tempo_index = []

for arquivo in os.listdir(f'./DensityMatrices/c_{cmod[0]}'):

    if rho_name in arquivo:
    
        tempo_index.append(int(arquivo.replace(rho_name, '').replace('.txt', '')))

tempo_index.sort()


for c in cmod:

    print(c)

    rhot_list = []

    ## reading rho(t)
    
    for t in tempo_index:

        rhot = np.loadtxt(f'./DensityMatrices/c_{c}/{rho_name}{t}.txt', dtype='complex', unpack=True)
    
        rhot_list.append(Qobj(rhot))


    ## Compute QFI

    tempo_QFI, QFI = Quantum_Fisher_Information(rhot_list, tempo_real, [[2],[2]])


    ## save QFI in file

    f = open(f'./FisherInformation/{modo}/QFI_{qubit}_c{c}.txt', 'w')

    for t in range(len(tempo_QFI)):

        f.write(f'{tempo_QFI[t]} {QFI[t]}\n')
        
        
    f.close()










