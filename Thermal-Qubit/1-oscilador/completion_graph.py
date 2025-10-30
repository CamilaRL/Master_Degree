import numpy as np
import matplotlib.pyplot as plt
import os


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
colors = iter(cmap(np.linspace(0.01, 1, tot)))

for i, curva in enumerate(curvas):

    tlist, completion = np.loadtxt(f'./ThermalKinematics_{modo}/completion_{curva}.txt', unpack=True)
    
    cor = next(colors)
    
    plt.plot(tlist, completion, color=cor, label=f'|c| {cmod[i]:.2e}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')

plt.legend(loc='best', bbox_to_anchor=(1., 0.25, 0.25, 0.6))
plt.tight_layout()
plt.show()

