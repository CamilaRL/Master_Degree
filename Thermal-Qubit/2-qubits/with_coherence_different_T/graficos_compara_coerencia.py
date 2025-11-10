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

cor=['blue', 'red']

cmod_list = []

for modo in modoList:
    cmod = np.loadtxt(f'./DensityMatrices/cmod_{modo}.txt', unpack=True, ndmin=1)

    cmod.sort()
    
    cmod_list.append(cmod)


for i, modo, in enumerate(modoList):
 
    for c in cmod_list[i]:
    
        tlist, completion = Read_Completion(c, modo)
    
        plt.plot(tlist, completion, color=cor[i])

plt.xlabel('Time')
plt.ylabel('Degree of completion')
plt.title(modo)
plt.text(0.5, 0.2, 'Red - Heating\nBlue - Cooling', bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1))
plt.xlim(left=0.01)
plt.xscale('log')
plt.show()


for i, modo, in enumerate(modoList):
 
    for c in cmod_list[i]:
    
        tlist, velocity = Read_Velocidade(c, modo)
    
        plt.plot(tlist, velocity, color=cor[i])

plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title(modo)
plt.text(0.5, 0.2, 'Red - Heating\nBlue - Cooling', bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1))
plt.xscale('log')
plt.xlim(left=0.01)
plt.show()

    
    
