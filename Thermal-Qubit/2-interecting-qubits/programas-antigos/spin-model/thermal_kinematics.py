import numpy as np
import matplotlib.pyplot as plt
import os



def Integracao(ydata, xdata):

    dt = xdata[1] - xdata[0]
    
    L = np.cumsum(ydata) * dt

    return L


def WriteOutput(tlist, velocity, position, completion, c, g, J, delta, modo):
    
    
    path = f'./ThermalKinematics_J{J}_delta{delta}/{modo}/'
    
    velocity_file = open(path + f'velocity_{c}_g{g}.txt', 'w')
    position_file = open(path + f'position_{c}_g{g}.txt', 'w')
    completion_file = open(path + f'completion_{c}_g{g}.txt', 'w')
    
    for t in range(len(tlist)):
    
        velocity_file.write(f'{tlist[t]} {velocity[t]}\n')
        position_file.write(f'{tlist[t]} {position[t]}\n')
        completion_file.write(f'{tlist[t]} {completion[t]}\n')
    
    velocity_file.close()
    position_file.close()
    completion_file.close()



### MAIN ###

qubit = 'q2'
modo = 'Cooling'
c = 'cmin'
title = 'Zero Initial Coherence'

J = 1
delta = -2

#os.mkdir(f'./ThermalKinematics_J{J}_delta{delta}')
#os.mkdir(f'./ThermalKinematics_J{J}_delta{delta}/{modo}')

## reading coherences
gList = np.loadtxt(f'./DensityMatrices_J{J}_delta{delta}/T_final_qubit.txt', unpack=True, usecols=(0))

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, len(gList))))

fig = plt.figure(figsize=(12,5))

for g in gList:
    
    print(g)
    
    curve_path = f'./FisherInformation_J{J}_delta{delta}/{modo}/QFI_{qubit}_{c}_g{g}.txt'
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    vlist = []
    
    for qfi in QFI:
        
        vlist.append(np.sqrt(qfi))
    
    cor = next(colors)
    
    Llist = Integracao(vlist, tlist)
    
    degree_completion = Llist/Llist[-1]
    
    plt.subplot(131)
    plt.scatter(tlist, vlist, color=cor, s=1, label=f'g = {g}')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    #plt.xscale('log')
    
    plt.subplot(132)
    plt.scatter(tlist, Llist, color=cor, s=1, label=f'g = {g}')
    plt.xlabel('Time')
    plt.ylabel('Position')
    #plt.xscale('log')

    plt.subplot(133)
    plt.plot(tlist, degree_completion, color=cor, label=f'g = {g}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')
    #plt.xscale('log')
    
    WriteOutput(tlist, vlist, Llist, degree_completion, c, g, J, delta, modo)

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.suptitle(modo+' - '+title)
plt.tight_layout()
plt.show()




