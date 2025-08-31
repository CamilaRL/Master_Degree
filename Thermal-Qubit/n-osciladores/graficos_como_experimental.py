import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(i, modo):

    path = f'./ThermalKinematics_{modo}/completion_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(i, modo):

    path = f'./ThermalKinematics_{modo}/velocity_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Mod_Cvalue(i, modo):

    path = f'./FisherInformation_{modo}/c_curva_{i}.txt'
    
    cvalue = np.loadtxt(path, unpack=True, dtype=complex, ndmin=1)
    
    cmod = abs(cvalue[0])
    
    return cmod


def Ordenamento(dados, curvas):

    dados, curvas_ordenadas = (list(t) for t in zip(*sorted(zip(dados, curvas))))
    
    return dados, curvas_ordenadas

### MAIN ###


modoList = ['Resfriar', 'Aquecer']

cList = []
curvasList = []

for modo in modoList:


    iList = []
    cvalue = []
    
    for arquivo in os.listdir(f'./FisherInformation_{modo}/'):

        if os.path.isfile(os.path.join(f'./FisherInformation_{modo}/', arquivo)) and 'c_curva_' in arquivo:
            
            curva_i = int(arquivo.replace('c_curva_', '').replace('.txt', ''))
            
            if curva_i > 0:
            
                iList.append(curva_i)

    for i in iList:
        
        cvalue.append(Read_Mod_Cvalue(i, modo))
    
    cvalue, iList = Ordenamento(cvalue, iList)
    
    curvasList.append(iList)
    cList.append(cvalue)


for i, curva in enumerate(curvasList[1]):
    
    
    print(f'Curva {curva} (aquecer) e {curvasList[0][i]} (resfriar)')
    print(f'c {modoList[0]}: {cList[0][i]}')
    print(f'c {modoList[1]}: {cList[1][i]}')
    
    
    tlist_r, completion_r = Read_Completion(curvasList[0][i], modoList[0])
    tlist_a, completion_a = Read_Completion(curva, modoList[1])
    
    plt.plot(tlist_r, completion_r, color='blue', label='Cooling')
    plt.plot(tlist_a, completion_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Degree of completion')
    #plt.title(f'Curva {curva} (aquecer) e {curvasList[0][i]} (resfriar)')
    plt.title(f'|c| = {cList[0][i]:.3f}')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.show()
    
    tlist_r, velocidade_r = Read_Velocidade(curvasList[0][i], modoList[0])
    tlist_a, velocidade_a = Read_Velocidade(curva, modoList[1])
    
    plt.plot(tlist_r, velocidade_r, color='blue', label='Cooling')
    plt.plot(tlist_a, velocidade_a, color='red', label='Heating')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    #plt.title(f'Curva {curva} (aquecer) e {curvasList[0][i]} (resfriar)')
    plt.title(f'|c| = {cList[0][i]:.3f}')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.ylim(top=23)
    plt.show()
    
    
    
    
    
