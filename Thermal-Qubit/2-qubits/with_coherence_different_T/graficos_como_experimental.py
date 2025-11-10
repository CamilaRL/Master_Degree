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


modoList = ['Cooling', 'Heating']

cList = []

cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True, ndmin=1)

cList.append(min(cmod))
cList.append(max(cmod))


for i in range(len(cList)):
    
    tlist_r, completion_r = Read_Completion(cList[i], modoList[0])
    tlist_a, completion_a = Read_Completion(cList[i], modoList[1])
    
    plt.plot(tlist_r, completion_r, color='blue', label='Cooling')
    plt.plot(tlist_a, completion_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Degree of completion')
    plt.title(f'|c| = {cList[i]:.2f}')
    plt.legend(loc='upper left')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.show()
    
    tlist_r, velocidade_r = Read_Velocidade(cList[i], modoList[0])
    tlist_a, velocidade_a = Read_Velocidade(cList[i], modoList[1])
    
    plt.plot(tlist_r, velocidade_r, color='blue', label='Cooling')
    plt.plot(tlist_a, velocidade_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.title(f'|c| = {cList[i]:.2f}')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.ylim(top=26)
    plt.show()
    
    
    
    
    
