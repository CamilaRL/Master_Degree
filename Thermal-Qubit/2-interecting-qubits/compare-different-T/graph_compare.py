import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(c, modo, g):

    path = f'./g-{g}/ThermalKinematics/{modo}/completion_{c}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(c, modo, g):

    path = f'./g-{g}/ThermalKinematics/{modo}/velocity_{c}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion




### MAIN ###

glist = [0, 0.8]

modoList = ['Cooling', 'Heating']

titulo = ['Zero Initial Coherence', 'Maximum Initial Coherence']

cmod0 = np.loadtxt(f'./g-{glist[0]}/DensityMatrices/cmod.txt', unpack=True, ndmin=1)
cmod1 = np.loadtxt(f'./g-{glist[1]}/DensityMatrices/cmod.txt', unpack=True, ndmin=1)

cList0 = [min(cmod0), max(cmod0)]
cList1 = [min(cmod1), max(cmod1)]


##### compare g

for i in range(len(cList0)):
    
    tlist_r0, completion_r0 = Read_Completion(cList0[i], modoList[0], glist[0])
    tlist_a0, completion_a0 = Read_Completion(cList0[i], modoList[1], glist[0])
    
    tlist_r1, completion_r1 = Read_Completion(cList1[i], modoList[0], glist[1])
    tlist_a1, completion_a1 = Read_Completion(cList1[i], modoList[1], glist[1])
    
    plt.plot(tlist_a0, completion_a0, color='red', linestyle='--', linewidth=2, label='g = 0 - Qubit 1')
    plt.plot(tlist_r0, completion_r0, color='blue', linestyle='--', linewidth=2, label='g = 0 - Qubit 2')
    
    plt.plot(tlist_a1, completion_a1, color='red', linewidth=2, label='g = 0.8 - Qubit 1')
    plt.plot(tlist_r1, completion_r1, color='blue', linewidth=2, label='g = 0.8 - Qubit 2')
    
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Degree of completion', fontsize=12)
    plt.title(titulo[i])
    plt.legend(loc='upper left', fontsize=12)
    plt.xscale('log')
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlim(left=0.1)
    plt.show()
    
    tlist_r0, velocidade_r0 = Read_Velocidade(cList0[i], modoList[0], glist[0])
    tlist_a0, velocidade_a0 = Read_Velocidade(cList0[i], modoList[1], glist[0])
    
    tlist_r1, velocidade_r1 = Read_Velocidade(cList1[i], modoList[0], glist[1])
    tlist_a1, velocidade_a1 = Read_Velocidade(cList1[i], modoList[1], glist[1])
    
    plt.plot(tlist_a0, velocidade_a0, color='red', linestyle='--', linewidth=2, label='g = 0 - Qubit 1')
    plt.plot(tlist_r0, velocidade_r0, color='blue', linestyle='--', linewidth=2, label='g = 0 - Qubit 2')
    
    plt.plot(tlist_a1, velocidade_a1, color='red', linewidth=2, label='g = 0.8 - Qubit 1')
    plt.plot(tlist_r1, velocidade_r1, color='blue', linewidth=2, label='g = 0.8 - Qubit 2')
    
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Velocity', fontsize=12)
    plt.title(titulo[i])
    plt.legend(loc='center right', fontsize=12)
    plt.xscale('log')
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlim(left=0.1)
    plt.ylim(top=0.8)
    plt.show()
    
    
    t_sigma_a0, sigma_a0 = np.loadtxt(f'./g-{glist[0]}/DensityMatrices/dSr_q1_c{cList0[i]}.txt', unpack=True)
    t_sigma_r0, sigma_r0 = np.loadtxt(f'./g-{glist[0]}/DensityMatrices/dSr_q2_c{cList0[i]}.txt', unpack=True)
    
    t_sigma_am, sigma_am = np.loadtxt(f'./g-{glist[1]}/DensityMatrices/dSr_q1_c{cList1[i]}.txt', unpack=True)
    t_sigma_rm, sigma_rm = np.loadtxt(f'./g-{glist[1]}/DensityMatrices/dSr_q2_c{cList1[i]}.txt', unpack=True)
    
    plt.plot(t_sigma_a0, sigma_a0, color='red', linestyle='--', linewidth=2, label='g = 0 - Qubit 1')
    plt.plot(t_sigma_r0, sigma_r0, color='blue', linestyle='--', linewidth=2, label='g = 0 - Qubit 2')
    
    plt.plot(t_sigma_am, sigma_am, color='red', linewidth=2, label='g = 0.8 - Qubit 1')
    plt.plot(t_sigma_rm, sigma_rm, color='blue', linewidth=2, label='g = 0.8 - Qubit 2')
    
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Local Entropy Production Rate', fontsize=12)
    plt.title(titulo[i])
    plt.legend(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xscale('log')
    plt.tight_layout()
    plt.show()
    
    

##### compare coherence

tlist_r0, completion_r0 = Read_Completion(cList0[0], modoList[0], glist[0])
tlist_a0, completion_a0 = Read_Completion(cList0[0], modoList[1], glist[0])

tlist_rmax, completion_rmax = Read_Completion(cList0[1], modoList[0], glist[0])
tlist_amax, completion_amax = Read_Completion(cList0[1], modoList[1], glist[0])

plt.plot(tlist_a0, completion_a0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cList0[0]:.3f} - Qubit 1')
plt.plot(tlist_r0, completion_r0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cList0[0]:.3f} - Qubit 2')

plt.plot(tlist_amax, completion_amax, color='orange', linestyle='--', linewidth=2, label=f'|c| = {cList0[1]:.3f} - Qubit 1')
plt.plot(tlist_rmax, completion_rmax, color='skyblue', linestyle='--', linewidth=2, label=f'|c| = {cList0[1]:.3f} - Qubit 2')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of completion', fontsize=12)
plt.title(f'g = 0')
plt.legend(loc='upper left', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)
plt.show()


tlist_r0, velocidade_r0 = Read_Velocidade(cList0[0], modoList[0], glist[0])
tlist_a0, velocidade_a0 = Read_Velocidade(cList0[0], modoList[1], glist[0])

tlist_rmax, velocidade_rmax = Read_Velocidade(cList0[1], modoList[0], glist[0])
tlist_amax, velocidade_amax = Read_Velocidade(cList0[1], modoList[1], glist[0])

plt.plot(tlist_a0, velocidade_a0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cList0[0]:.3f} - Qubit 1')
plt.plot(tlist_r0, velocidade_r0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cList0[0]:.3f} - Qubit 2')

plt.plot(tlist_amax, velocidade_amax, color='orange', linestyle='--', linewidth=2, label=f'|c| = {cList0[1]:.3f} - Qubit 1')
plt.plot(tlist_rmax, velocidade_rmax, color='skyblue', linestyle='--', linewidth=2, label=f'|c| = {cList0[1]:.3f} - Qubit 2')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Velocity', fontsize=12)
plt.title(f'g = 0')
plt.legend(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)
plt.show()




tlist_r0, completion_r0 = Read_Completion(cList1[0], modoList[0], glist[1])
tlist_a0, completion_a0 = Read_Completion(cList1[0], modoList[1], glist[1])

tlist_rmax, completion_rmax = Read_Completion(cList1[1], modoList[0], glist[1])
tlist_amax, completion_amax = Read_Completion(cList1[1], modoList[1], glist[1])

plt.plot(tlist_a0, completion_a0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cList1[0]:.3f} - Qubit 1')
plt.plot(tlist_r0, completion_r0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cList1[0]:.3f} - Qubit 2')

plt.plot(tlist_amax, completion_amax, color='orange', linestyle='-', linewidth=2, label=f'|c| = {cList1[1]:.3f} - Qubit 1')
plt.plot(tlist_rmax, completion_rmax, color='skyblue', linestyle='-', linewidth=2, label=f'|c| = {cList1[1]:.3f} - Qubit 2')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of completion', fontsize=12)
plt.title(f'g = 0.8')
plt.legend(loc='upper left', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)
plt.show()



tlist_r0, velocidade_r0 = Read_Velocidade(cList1[0], modoList[0], glist[1])
tlist_a0, velocidade_a0 = Read_Velocidade(cList1[0], modoList[1], glist[1])

tlist_rmax, velocidade_rmax = Read_Velocidade(cList1[1], modoList[0], glist[1])
tlist_amax, velocidade_amax = Read_Velocidade(cList1[1], modoList[1], glist[1])

plt.plot(tlist_a0, velocidade_a0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cList1[0]:.3f} - Qubit 1')
plt.plot(tlist_r0, velocidade_r0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cList1[0]:.3f} - Qubit 2')

plt.plot(tlist_amax, velocidade_amax, color='orange', linestyle='-', linewidth=2, label=f'|c| = {cList1[1]:.3f} - Qubit 1')
plt.plot(tlist_rmax, velocidade_rmax, color='skyblue', linestyle='-', linewidth=2, label=f'|c| = {cList1[1]:.3f} - Qubit 2')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Velocity', fontsize=12)
plt.title(f'g = 0.8')
plt.legend(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)
plt.show()

