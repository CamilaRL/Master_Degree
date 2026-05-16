import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(mol, c, modo, g):

    path = f'./{mol}-model/ThermalKinematics/{modo}/completion_{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(mol, c, modo, g):

    path = f'./{mol}-model/ThermalKinematics/{modo}/velocity_{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion




### MAIN ###

modelos = ['XXX', 'XXZ', 'XYZ', 'XX', 'XY', 'TFI']

modoList = ['Cooling', 'Heating']



##### compare g

plt.figure(figsize=(10,5))

for mol in modelos:
    
    glist = np.loadtxt(f'./{mol}-model/DensityMatrices/correlations.txt', unpack=True, usecols=(0))
    
    J = glist[1]
    
    tlist_r, completion_r = Read_Completion(mol, 'cmax', modoList[0], J)
    tlist_a, completion_a = Read_Completion(mol, 'cmax', modoList[1], J)
    
    plt.subplot(121)
    
    plt.plot(tlist_a, completion_a, linewidth=2, label=f'{mol}')
    
    plt.ylabel('Degree of Completion', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.title('Heating')
    plt.legend(loc='center right', fontsize=12)
    
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
    plt.subplot(122)
    plt.plot(tlist_r, completion_r, linewidth=2, label=f'{mol}')
    
    plt.xlabel('Time', fontsize=12)
    plt.title('Cooling')
    plt.legend(loc='center right', fontsize=12)
    
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)

plt.show()


for mol in modelos:
    
    glist = np.loadtxt(f'./{mol}-model/DensityMatrices/correlations.txt', unpack=True, usecols=(0))
    
    J = glist[1]
    
    tlist_r, completion_r = Read_Completion(mol, 'cmax', modoList[0], J)
    tlist_a, completion_a = Read_Completion(mol, 'cmax', modoList[1], J)
    
    dif_completion = []

    for t in range(len(tlist_a)):
        
        dif_completion.append(completion_a[t] - completion_r[t])
    
    
    plt.plot(tlist_a, dif_completion, linewidth=2, label=f'{mol}')
    
plt.ylabel('Degree of Completion Difference', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.legend(loc='center right', fontsize=12)

plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.show()
    