import numpy as np
import matplotlib.pyplot as plt
import os



def Integracao(ydata, xdata):

    dt = xdata[1] - xdata[0]
    
    L = np.cumsum(ydata) * dt

    return L


def WriteOutput(tlist, velocity, position, completion, i):
    
    path = f'./ThermalKinematics/'
    
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

modo = 'Aquecer'
cmod = []
curvas = []

for arquivo in os.listdir(f'./FisherInformation_{modo}/'):

    path_to_file = os.path.join(f'./FisherInformation_{modo}/', arquivo)
    
    if os.path.isfile(path_to_file) and 'c_curva_' in arquivo:
    
        curvas.append(int(arquivo.replace('c_curva_', '').replace('.txt', '')))
        
        cList = np.loadtxt(path_to_file, unpack=True, ndmin=1, dtype='complex')
        
        cmod.append(abs(cList[0]))

cmod, curvas = (list(t) for t in zip(*sorted(zip(cmod, curvas))))


tot = len(cmod)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, tot)))
fig = plt.figure(figsize=(12,5))

for i, curva in enumerate(curvas):
    
    curve_path = f'./FisherInformation_{modo}/curva_{curva}.txt'
    
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    vlist = []
    
    for qfi in QFI:
        
        vlist.append(np.sqrt(qfi))
    
    cor = next(colors)
    
    Llist = Integracao(vlist, tlist)
    
    degree_completion = Llist/Llist[-1]
    
    plt.subplot(131)
    plt.scatter(tlist, vlist, color=cor, s=1, label=f'|c| {cmod[i]:.2e}')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    
    plt.subplot(132)
    plt.scatter(tlist, Llist, color=cor, s=1, label=f'|c| {cmod[i]:.2e}')
    plt.xlabel('Time')
    plt.ylabel('Position')

    plt.subplot(133)
    plt.plot(tlist, degree_completion, color=cor, label=f'|c| {cmod[i]:.2e}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')
    
    WriteOutput(tlist, vlist, Llist, degree_completion, i)

plt.legend(loc='best', bbox_to_anchor=(1., 0.25, 0.25, 0.6))
plt.tight_layout()
plt.show()



