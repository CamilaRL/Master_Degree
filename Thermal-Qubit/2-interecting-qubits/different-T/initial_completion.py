import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(c, modo):

    path = f'./ThermalKinematics/{modo}/completion_{c}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(c, modo):

    path = f'./ThermalKinematics/{modo}/velocity_{c}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion
    
    
    
### MAIN ###


modo = 'Heating'
cor = 'red'

cList = []

cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True, ndmin=1)
    
cmod.sort()

initial_completion = []
    
for c in cmod:
    
    tlist, completion = Read_Completion(c, modo)
        
    initial_completion.append(completion[0])

plt.scatter(cmod, initial_completion, color=cor, label=modo)
plt.plot(cmod, initial_completion, color=cor)

plt.legend()
plt.ylabel('Initial Completion')
plt.xlabel('|c|')
plt.tight_layout()
plt.show()


