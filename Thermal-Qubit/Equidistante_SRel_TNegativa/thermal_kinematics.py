import numpy as np
import matplotlib.pyplot as plt
import os

def Integracao(ydata, xdata):

    dt = xdata[1] - xdata[0]
    
    L = np.cumsum(ydata) * dt

    return L


def WriteOutput(tlist, velocity, position, completion, i, modo, dSr):
    
    
    path = f'./ThermalKinematics_{modo}_{dSr}/'
    
    velocity_file = open(path + f'velocity_{i}.txt', 'w')
    position_file = open(path + f'position_{i}.txt', 'w')
    completion_file = open(path + f'completion_{i}.txt', 'w')
    
    for t in range(len(tlist)):
    
        velocity_file.write(f'{tlist[t]} {velocity[t]}\n')
        position_file.write(f'{tlist[t]} {position[t]}\n')
        completion_file.write(f'{tlist[t]} {completion[t]}\n')
    
    velocity_file.close()
    position_file.close()
    completion_file.close()



### MAIN ###

modo = 'Cooling'
dSr = 0.1

curvas = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)

os.mkdir(f'./ThermalKinematics_{modo}_{dSr}')

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, len(curvas))))

fig = plt.figure(figsize=(12,5))


for i, j in enumerate(curvas):
    
    curva = int(j)
    
    curve_path = f'./FisherInformation_{modo}_{dSr}/QFI_curve_{curva}.txt'
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    vlist = []
    
    for qfi in QFI:
        
        vlist.append(np.sqrt(qfi))
    
    cor = next(colors)
    
    Llist = Integracao(vlist, tlist)
    
    degree_completion = Llist/Llist[-1]
    
    plt.subplot(131)
    plt.scatter(tlist, vlist, color=cor, s=1, label=f'|c| = {cmod[i]:.3f}')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    
    plt.subplot(132)
    plt.scatter(tlist, Llist, color=cor, s=1, label=f'|c| = {cmod[i]:.3f}')
    plt.xlabel('Time')
    plt.ylabel('Position')

    plt.subplot(133)
    plt.plot(tlist, degree_completion, color=cor, label=f'|c| = {cmod[i]:.3f}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')
    
    WriteOutput(tlist, vlist, Llist, degree_completion, curva, modo, dSr)

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.suptitle(modo)
plt.tight_layout()
plt.show()



