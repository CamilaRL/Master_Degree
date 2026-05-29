import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(c, modo, g):

    path = f'./ThermalKinematics/{modo}/completion_c{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(c, modo, g):

    path = f'./ThermalKinematics/{modo}/velocity_c{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion




### MAIN ###

modoList = ['Cooling', 'Heating']

titulo = ['Zero Initial Coherence', 'Maximum Initial Coherence']

cList = ['min', 'max']

glist, cmin, cmax = np.loadtxt(f'DensityMatrices/correlations.txt', unpack=True)


##### compare g
fig = plt.figure(figsize=(10,5))

for i in range(len(cList)):
    
    num = 121 + i
    
    plt.subplot(num)
    
    tlist_r0, completion_r0 = Read_Completion(cList[i], modoList[0], glist[0])
    tlist_a0, completion_a0 = Read_Completion(cList[i], modoList[1], glist[0])
    
    tlist_r1, completion_r1 = Read_Completion(cList[i], modoList[0], glist[1])
    tlist_a1, completion_a1 = Read_Completion(cList[i], modoList[1], glist[1])
    
    plt.plot(tlist_a0, completion_a0, color='red', linestyle='--', linewidth=2, label='J = 0 - Heating')
    plt.plot(tlist_r0, completion_r0, color='blue', linestyle='--', linewidth=2, label='J = 0 - Cooling')
    
    plt.plot(tlist_a1, completion_a1, color='red', linewidth=2, label='J = 0.8 - Heating')
    plt.plot(tlist_r1, completion_r1, color='blue', linewidth=2, label='J = 0.8 - Cooling')
    
    plt.xlabel('Time', fontsize=12)
    plt.title(titulo[i], fontsize=14)
    plt.legend(loc='center right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
    if i == 0:
        plt.ylabel('Degree of Completion', fontsize=12)

plt.tight_layout()
plt.show()

fig = plt.figure(figsize=(10,5))

for i in range(len(cList)):
    
    num = 121 + i
    
    plt.subplot(num)
    
    tlist_r0, velocidade_r0 = Read_Velocidade(cList[i], modoList[0], glist[0])
    tlist_a0, velocidade_a0 = Read_Velocidade(cList[i], modoList[1], glist[0])
    
    tlist_r1, velocidade_r1 = Read_Velocidade(cList[i], modoList[0], glist[1])
    tlist_a1, velocidade_a1 = Read_Velocidade(cList[i], modoList[1], glist[1])
    
    plt.plot(tlist_a0, velocidade_a0, color='red', linestyle='--', linewidth=2, label='J = 0 - Heating')
    plt.plot(tlist_r0, velocidade_r0, color='blue', linestyle='--', linewidth=2, label='J = 0 - Cooling')
    
    plt.plot(tlist_a1, velocidade_a1, color='red', linewidth=2, label='J = 0.8 - Heating')
    plt.plot(tlist_r1, velocidade_r1, color='blue', linewidth=2, label='J = 0.8 - Cooling')
    
    plt.xlabel('Time', fontsize=12)
    plt.title(titulo[i], fontsize=14)
    plt.legend(loc='center right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlim(right=18)
    
    if i == 0:
        plt.ylabel('Velocity', fontsize=12)

plt.tight_layout()
plt.show()
    
fig = plt.figure(figsize=(10,5))

for i in range(len(cList)):
    
    num = 121 + i
    
    plt.subplot(num)
    
    t_sigma_a0, sigma_a0 = np.loadtxt(f'./Thermodynamics/dSr_q1_c{cList[i]}_g{glist[0]}.txt', unpack=True)
    t_sigma_r0, sigma_r0 = np.loadtxt(f'./Thermodynamics/dSr_q2_c{cList[i]}_g{glist[0]}.txt', unpack=True)
    
    t_sigma_am, sigma_am = np.loadtxt(f'./Thermodynamics/dSr_q1_c{cList[i]}_g{glist[1]}.txt', unpack=True)
    t_sigma_rm, sigma_rm = np.loadtxt(f'./Thermodynamics/dSr_q2_c{cList[i]}_g{glist[1]}.txt', unpack=True)
    
    plt.plot(t_sigma_a0, sigma_a0, color='red', linestyle='--', linewidth=2, label='J = 0 - Heating')
    plt.plot(t_sigma_r0, sigma_r0, color='blue', linestyle='--', linewidth=2, label='J = 0 - Cooling')
    
    plt.plot(t_sigma_am, sigma_am, color='red', linewidth=2, label='J = 0.8 - Heating')
    plt.plot(t_sigma_rm, sigma_rm, color='blue', linewidth=2, label='J = 0.8 - Cooling')
    
    plt.xlabel('Time', fontsize=12)
    plt.title(titulo[i], fontsize=14)
    plt.legend(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlim(right=15)
    
    if i == 0:
        plt.ylabel('Local Entropy Production Rate', fontsize=12)
    
plt.tight_layout()
plt.show()
    


##### compare coherence
for i, g in enumerate(glist):
    
    fig = plt.figure(figsize=(10,5))
    
    plt.subplot(121)
    
    tlist_r0, completion_r0 = Read_Completion(cList[0], modoList[0], g)
    tlist_a0, completion_a0 = Read_Completion(cList[0], modoList[1], g)

    tlist_rmax, completion_rmax = Read_Completion(cList[1], modoList[0], g)
    tlist_amax, completion_amax = Read_Completion(cList[1], modoList[1], g)

    plt.plot(tlist_a0, completion_a0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cmin[i]:.3f} - Heating')
    plt.plot(tlist_r0, completion_r0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cmin[i]:.3f} - Cooling')

    plt.plot(tlist_amax, completion_amax, color='red', linestyle='--', linewidth=2, label=f'|c| = {cmax[i]:.3f} - Heating')
    plt.plot(tlist_rmax, completion_rmax, color='blue', linestyle='--', linewidth=2, label=f'|c| = {cmax[i]:.3f} - Cooling')

    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Degree of Completion', fontsize=12)
    plt.legend(loc='lower right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xscale('log')
    plt.xlim(left=0.1)

    plt.subplot(122)
    
    tlist_r0, velocidade_r0 = Read_Velocidade(cList[0], modoList[0], g)
    tlist_a0, velocidade_a0 = Read_Velocidade(cList[0], modoList[1], g)

    tlist_rmax, velocidade_rmax = Read_Velocidade(cList[1], modoList[0], g)
    tlist_amax, velocidade_amax = Read_Velocidade(cList[1], modoList[1], g)

    plt.plot(tlist_a0, velocidade_a0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cmin[i]:.3f} - Heating')
    plt.plot(tlist_r0, velocidade_r0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cmin[i]:.3f} - Cooling')

    plt.plot(tlist_amax, velocidade_amax, color='red', linestyle='--', linewidth=2, label=f'|c| = {cmax[i]:.3f} - Heating')
    plt.plot(tlist_rmax, velocidade_rmax, color='blue', linestyle='--', linewidth=2, label=f'|c| = {cmax[i]:.3f} - Cooling')

    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Velocity', fontsize=12)
    plt.legend(loc='upper right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xscale('log')
    plt.xlim(left=0.1)
    
    plt.suptitle(f'J = {g}', fontsize=14)
    plt.tight_layout()
    plt.show()

