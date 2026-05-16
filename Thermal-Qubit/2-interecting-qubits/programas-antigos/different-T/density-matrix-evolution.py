import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
import os


def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))


def Interseccao_Temperatura(p):
    
    find_T = lambda T: pFunc(T, w0=2) - p
        
    T = brentq(find_T, 0.01, 50)  # busca em todo o range de T
    
    return T


####### MAIN #######

modo = 'Heating'
rhoname = 'rhof_q1_t'

cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)

cmod_extremes = [min(cmod), max(cmod)]

tlist = np.arange(0.005, 20, 0.005)


for c in cmod_extremes:

    Tt = []
    coherence = []
    timestep = []

    path = f'./DensityMatrices/c_{c}/'

    for arquivo in os.listdir(path):
    
        if rhoname in arquivo:
           
           timestep.append(int(arquivo.replace(rhoname, '').replace('.txt', '')))
            
    timestep.sort()
    
    for t in timestep:
    
        rhot = np.loadtxt(path + rhoname + f'{t}.txt', unpack=True, dtype='complex')
            
        #print(Qobj(rhot).tr(), Qobj(rhot).isherm)
        
        T = Interseccao_Temperatura(rhot[0][0].real)
            
        Tt.append(T)
            
        #coherence.append( abs(rhot[0][1]) )
                        
    plt.plot(tlist, Tt, label=f'|c| = {c:.6f}')
        
plt.xscale('log')
plt.ylabel('Temperature')
plt.xlabel('Time')
plt.title(modo)
plt.legend()
plt.show()

########################## TOTAL COHERENCE #############################

rhoname = 'rhof_t'

for c in cmod_extremes:

    coherence = []
    timestep = []

    path = f'./DensityMatrices/c_{c}/'

    for arquivo in os.listdir(path):
    
        if rhoname in arquivo:
           
           timestep.append(int(arquivo.replace(rhoname, '').replace('.txt', '')))
            
    timestep.sort()
    
    for t in timestep:

        rhot = np.loadtxt(path + rhoname + f'{t}.txt', unpack=True, dtype='complex')
        
        coherence.append(abs(rhot[1][2]))
        
    
    plt.plot(tlist, coherence, label=f'|c| = {c:.6f}')

plt.xscale('log')
plt.ylabel('Total Coherence')
plt.xlabel('Time')
plt.title(modo)
plt.legend()
plt.show()
