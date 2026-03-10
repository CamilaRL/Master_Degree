import numpy as np
import matplotlib.pyplot as plt
import os



def Integracao(ydata, xdata):

    dt = xdata[1] - xdata[0]
    
    L = np.cumsum(ydata) * dt

    return L


def WriteOutput(tlist, velocity, position, completion, c, g, modo):
    
    
    path = f'./ThermalKinematics/{modo}/'
    
    velocity_file = open(path + f'velocity_c{c}_g{g}.txt', 'w')
    position_file = open(path + f'position_c{c}_g{g}.txt', 'w')
    completion_file = open(path + f'completion_c{c}_g{g}.txt', 'w')
    
    for t in range(len(tlist)):
    
        velocity_file.write(f'{tlist[t]} {velocity[t]}\n')
        position_file.write(f'{tlist[t]} {position[t]}\n')
        completion_file.write(f'{tlist[t]} {completion[t]}\n')
    
    velocity_file.close()
    position_file.close()
    completion_file.close()



### MAIN ###

qubit = 'q1'
modo = 'Heating'

os.mkdir(f'./ThermalKinematics/{modo}')


gList = [0.0, 0.8]

cList = ['min', 'max']

fig = plt.figure(figsize=(12,5))

for c in cList:

    for g in gList:
    
        tlist, QFI = np.loadtxt(f'./FisherInformation/{modo}/QFI_{qubit}_c{c}_g{g}.txt', unpack=True)
    
        vlist = []
    
        for qfi in QFI:
        
            vlist.append(np.sqrt(qfi))
        
        Llist = Integracao(vlist, tlist)
        
        degree_completion = Llist/Llist[-1]
        
        plt.subplot(131)
        plt.plot(tlist, vlist, label=f'c{c} | g = {g}')
        plt.xlabel('Time')
        plt.ylabel('Velocity')
        plt.xscale('log')
        
        plt.subplot(132)
        plt.plot(tlist, Llist, label=f'c{c} | g = {g}')
        plt.xlabel('Time')
        plt.ylabel('Position')
        plt.xscale('log')

        plt.subplot(133)
        plt.plot(tlist, degree_completion, label=f'c{c} | g = {g}')
        plt.xlabel('Time')
        plt.ylabel('Degree of Completion')
        plt.xscale('log')
        
        WriteOutput(tlist, vlist, Llist, degree_completion, c, g, modo)


plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.suptitle(modo)
plt.tight_layout()
plt.show()
