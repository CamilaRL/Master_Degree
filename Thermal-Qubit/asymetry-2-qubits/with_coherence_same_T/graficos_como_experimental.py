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

for modo in modoList:

    cmod = np.loadtxt(f'./DensityMatrices/cmod_{modo}.txt', unpack=True, ndmin=1)
    
    cList.append([min(cmod), max(cmod)])


for i in range(len(cList[0])):
    
    print(f'c {modoList[0]}: {cList[0][i]}')
    print(f'c {modoList[1]}: {cList[1][i]}')
    
    
    tlist_r, completion_r = Read_Completion(cList[0][i], modoList[0])
    tlist_a, completion_a = Read_Completion(cList[1][i], modoList[1])
    
    plt.plot(tlist_r, completion_r, color='blue', label='Cooling')
    plt.plot(tlist_a, completion_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Degree of completion')
    plt.title(r'|$c_{cooling}$| = ' + f'{cList[0][i]:.5f} and ' r'|$c_{heating}$| = ' + f'{cList[1][i]:.5f}')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.show()
    
    tlist_r, velocidade_r = Read_Velocidade(cList[0][i], modoList[0])
    tlist_a, velocidade_a = Read_Velocidade(cList[1][i], modoList[1])
    
    plt.plot(tlist_r, velocidade_r, color='blue', label='Cooling')
    plt.plot(tlist_a, velocidade_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.title(r'|$c_{cooling}$| = ' + f'{cList[0][i]:.5f} and ' r'|$c_{heating}$| = ' + f'{cList[1][i]:.5f}')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.ylim(top=26)
    plt.show()
    
    
    
    
    
