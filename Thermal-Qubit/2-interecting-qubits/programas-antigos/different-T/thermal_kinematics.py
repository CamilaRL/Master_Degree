import numpy as np
import matplotlib.pyplot as plt
import os



def Integracao(ydata, xdata):

    dt = xdata[1] - xdata[0]
    
    L = np.cumsum(ydata) * dt

    return L


def WriteOutput(tlist, velocity, position, completion, c, modo):
    
    
    path = f'./ThermalKinematics/{modo}/'
    
    velocity_file = open(path + f'velocity_{c}.txt', 'w')
    position_file = open(path + f'position_{c}.txt', 'w')
    completion_file = open(path + f'completion_{c}.txt', 'w')
    
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


## reading coherences
cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, len(cmod))))

fig = plt.figure(figsize=(12,5))

for i, c in enumerate(cmod):
    
    print(i)
    
    curve_path = f'./FisherInformation/{modo}/QFI_{qubit}_c{c}.txt'
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    vlist = []
    
    for qfi in QFI:
        
        vlist.append(np.sqrt(qfi))
    
    cor = next(colors)
    
    Llist = Integracao(vlist, tlist)
    
    degree_completion = Llist/Llist[-1]
    
    plt.subplot(131)
    plt.scatter(tlist, vlist, color=cor, s=1, label=f'|c| = {c}')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.xscale('log')
    
    plt.subplot(132)
    plt.scatter(tlist, Llist, color=cor, s=1, label=f'|c| = {c}')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.xscale('log')

    plt.subplot(133)
    plt.plot(tlist, degree_completion, color=cor, label=f'|c| = {c}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')
    plt.xscale('log')
    
    WriteOutput(tlist, vlist, Llist, degree_completion, c, modo)

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.suptitle(modo)
plt.tight_layout()
plt.show()
