import numpy as np
import matplotlib.pyplot as plt
from qutip import *

def Derivada_Rho(rho_list, tlist):
    
    drho_list = [np.zeros((2, 2), dtype=complex) for k in range(len(tlist)-1)]
    
    for i in range(2):
        for j in range(2):
            
            for t in range(len(tlist)-1):
    
                rho_i = rho_list[t][i][j]
                rho_f = rho_list[t+1][i][j]
                temp_i = tlist[t]
                temp_f = tlist[t+1]

                drho_list[t][i][j] = (rho_f - rho_i)/(temp_f - temp_i)
                
    return drho_list
    

## parameters

w0 = 1
g = 0


## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2



## Heat

cmod = np.loadtxt(f'./g-{g}/DensityMatrices/cmod.txt', unpack=True)

cmod_extremes = [min(cmod), max(cmod)]

tlist = np.arange(0, 20, 0.01)
timestep = np.arange(0, len(tlist), 1)

fig = plt.figure(figsize=(10,5))

cor = ['red', 'orange', 'blue', 'skyblue']
linha = ['-', '--', '-', '--']

for k, c in enumerate(cmod_extremes):
    
    print(c)
    
    path = f'./g-{g}/DensityMatrices/'
    
    rho1_list = []
    rho2_list = []
    
    for t in timestep:
        
        rhot1 = np.loadtxt(path + f'c_{c}/rhof_q1_t{t}.txt', unpack=True, dtype='complex')
        rhot2 = np.loadtxt(path + f'c_{c}/rhof_q2_t{t}.txt', unpack=True, dtype='complex')
    
        rho1_list.append(rhot1)
        rho2_list.append(rhot2)
    
    drho1 = Derivada_Rho(rho1_list, tlist)
    drho2 = Derivada_Rho(rho2_list, tlist)
    
    Q1 = []
    Q2 = []
    
    for t in range(len(timestep)-1):
        
        Q1.append( (Qobj(drho1[t]) * H_q1).tr() )
        Q2.append( (Qobj(drho2[t]) * H_q2).tr() )
        
        
    f1 = open(path + f'Q1_c{c}.txt', 'w')
    f2 = open(path + f'Q2_c{c}.txt', 'w')

    for t in range(len(tlist)-1):
        
        f1.write(f'{tlist[t]} {Q1[t]}\n')
        f2.write(f'{tlist[t]} {Q2[t]}\n')
        
    f1.close()
    f2.close()
        
    plt.subplot(121)
    plt.plot(tlist[:-1], Q1, color=cor[k], linestyle=linha[k], label=f'|c| = {c:.3f}')
    
    plt.subplot(122)
    plt.plot(tlist[:-1], Q2, color=cor[k+2], linestyle=linha[k+2], label=f'|c| = {c:.3f}')

plt.subplot(121)
plt.ylabel('Heat Flow')
plt.xlabel('Time')
plt.title('Qubit 1')
plt.legend()

plt.subplot(122)
plt.ylabel('Heat Flow')
plt.xlabel('Time')
plt.title('Qubit 2')

plt.suptitle(f'g = {g}')
plt.legend()
plt.tight_layout()
plt.show()