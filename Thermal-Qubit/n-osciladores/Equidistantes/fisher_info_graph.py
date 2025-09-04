import numpy as np
import matplotlib.pyplot as plt
import os


def Read_Mod_Cvalue(i, dT, modo):

    path = f'./dT-{dT}/FisherInformation_{modo}/c_curva_{i}.txt'
    
    cvalue = np.loadtxt(path, unpack=True, dtype=complex, ndmin=1)
    
    cmod = abs(cvalue[0])
    
    return cmod


def Ordenamento(dados, curvas):

    dados, curvas_ordenadas = (list(t) for t in zip(*sorted(zip(dados, curvas))))
    
    return dados, curvas_ordenadas
    
    

### MAIN ###

modo = 'Resfriar'
dT = 1

iList = []
cvalue = []
    
for arquivo in os.listdir(f'./dT-{dT}/FisherInformation_{modo}/'):

    if os.path.isfile(os.path.join(f'./dT-{dT}/FisherInformation_{modo}/', arquivo)) and 'c_curva_' in arquivo:
            
        curva_i = int(arquivo.replace('c_curva_', '').replace('.txt', ''))
            
        if curva_i > 0:
            
            iList.append(curva_i)

for i in iList:
        
    cvalue.append(Read_Mod_Cvalue(i, dT, modo))
    
cvalue, iList = Ordenamento(cvalue, iList)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(iList))))

curvas = []

for i in range(len(iList)):

    tlist, QFI = np.loadtxt(f'./dT-{dT}/FisherInformation_{modo}/curva_{iList[i]}.txt', unpack=True)
    
    cor = next(colors)
    
    plt.scatter(tlist, QFI, color=cor, s=10, label=f'|c| = {abs(cvalue[i]):.3f}')
    plt.xlabel('Time')
    plt.ylabel('Quantum Fisher Information')

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.title(r'Cooling ($\Delta$T = ' + f'{dT})')
plt.xscale('log')
plt.tight_layout()
plt.show()
