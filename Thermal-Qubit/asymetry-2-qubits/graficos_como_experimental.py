import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(i, modo):

    path = f'./ThermalKinematics_{modo}/completion_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(i, modo):

    path = f'./ThermalKinematics_{modo}/velocity_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion




### MAIN ###

Sr = 0.1
modoList = ['Cooling', 'Heating']

cList = []
curvasList = []

'''for modo in modoList:

    curvas = np.loadtxt(f'./FisherInformation_{modo}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
    cmod = np.loadtxt(f'./FisherInformation_{modo}/cmod.txt', unpack=True, usecols=(1), ndmin=1)
    
    curvasList.append(curvas)
    cList.append(cmod)'''


for i, curva in enumerate([0]):
    
    curva_r = 0
    curva_a = 0
    
    '''print(f'Curva {curva} (aquecer) e {curvasList[0][i]} (resfriar)')
    print(f'c {modoList[0]}: {cList[0][i]}')
    print(f'c {modoList[1]}: {cList[1][i]}')'''
    
    
    tlist_r, completion_r = Read_Completion(curva_r, modoList[0])
    tlist_a, completion_a = Read_Completion(curva_a, modoList[1])
    
    plt.plot(tlist_r, completion_r, color='blue', label='Cooling')
    plt.plot(tlist_a, completion_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Degree of completion')
    #plt.title(f'Curva {curva} (aquecer) e {curvasList[0][i]} (resfriar)')
    plt.title(f'|c| = 0')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.show()
    
    tlist_r, velocidade_r = Read_Velocidade(curva_r, modoList[0])
    tlist_a, velocidade_a = Read_Velocidade(curva_a, modoList[1])
    
    plt.plot(tlist_r, velocidade_r, color='blue', label='Cooling')
    plt.plot(tlist_a, velocidade_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    #plt.title(f'Curva {curva} (aquecer) e {curvasList[0][i]} (resfriar)')
    plt.title(f'|c| = 0')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.ylim(top=26)
    plt.show()
    
    
    
    
    
