import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(c, modo, g, J, delta):

    path = f'./ThermalKinematics_J{J}_delta{delta}/{modo}/completion_{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(c, modo, g, J, delta):

    path = f'./ThermalKinematics_J{J}_delta{delta}/{modo}/velocity_{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion




### MAIN ###

J = 1
delta = -2

gList = np.loadtxt(f'./DensityMatrices_J{J}_delta{delta}/T_final_qubit.txt', unpack=True, usecols=(0))

modoList = ['Cooling', 'Heating']

cName = ['Zero Initial Coherence', 'Maximum Initial Coherence']

cList = ['cmin', 'cmax']


##### compare g

for g in gList:
    
    tlist_r_min, completion_r_min = Read_Completion(cList[0], modoList[0], g, J, delta)
    tlist_a_min, completion_a_min = Read_Completion(cList[0], modoList[1], g, J, delta)
    
    tlist_r_max, completion_r_max = Read_Completion(cList[1], modoList[0], g, J, delta)
    tlist_a_max, completion_a_max = Read_Completion(cList[1], modoList[1], g, J, delta)
    
    plt.plot(tlist_a_min, completion_a_min, color='red', linewidth=2, label=f'{cName[0]} - {modoList[1]}')
    plt.plot(tlist_r_min, completion_r_min, color='blue', linewidth=2, label=f'{cName[0]} - {modoList[0]}')
    
    plt.plot(tlist_a_max, completion_a_max, color='red', linestyle='--', linewidth=2, label=f'{cName[1]} - {modoList[1]}')
    plt.plot(tlist_r_max, completion_r_max, color='blue', linestyle='--', linewidth=2, label=f'{cName[1]} - {modoList[0]}')
    
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Degree of completion', fontsize=12)
    plt.title(f'g = {g}', fontsize=14)
    plt.legend(loc='upper left', fontsize=12)
    plt.xscale('log')
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    #plt.xlim(left=0.1)
    plt.tight_layout()
    plt.show()
    
    tlist_r_min, velocidade_r_min = Read_Velocidade(cList[0], modoList[0], g, J, delta)
    tlist_a_min, velocidade_a_min = Read_Velocidade(cList[0], modoList[1], g, J, delta)
    
    tlist_r_max, velocidade_r_max = Read_Velocidade(cList[1], modoList[0], g, J, delta)
    tlist_a_max, velocidade_a_max = Read_Velocidade(cList[1], modoList[1], g, J, delta)
    
    plt.plot(tlist_a_min, velocidade_a_min, color='red', linewidth=2, label=f'{cName[0]} - {modoList[1]}')
    plt.plot(tlist_r_min, velocidade_r_min, color='blue', linewidth=2, label=f'{cName[0]} - {modoList[0]}')
    
    plt.plot(tlist_a_max, velocidade_a_max, color='red', linestyle='--', linewidth=2, label=f'{cName[1]} - {modoList[1]}')
    plt.plot(tlist_r_max, velocidade_r_max, color='blue', linestyle='--', linewidth=2, label=f'{cName[1]} - {modoList[0]}')
    
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Velocity', fontsize=12)
    plt.title(f'g = {g}', fontsize=14)
    plt.legend(loc='best', fontsize=12)
    #plt.legend(loc='best', fontsize=12, bbox_to_anchor=(1., 0.5))
    plt.xscale('log')
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    #plt.xlim(left=0.1)
    plt.tight_layout()
    plt.show()


