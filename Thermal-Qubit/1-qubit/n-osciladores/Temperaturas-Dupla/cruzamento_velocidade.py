import numpy as np
import matplotlib.pyplot as plt
import os
import math

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
tlist_inter = []
va_inter = []  
vr_inter = []
Sa_inter = []
Sr_inter = []
Sa_max = []
Sr_max = []

for i, curva in enumerate(curvasList[1]):
    
    tlist, S_r = np.loadtxt(f'./Entropy_{modoList[0]}/entropy-{cList[0][i]:.3f}.txt', unpack=True)
    tlist, S_a = np.loadtxt(f'./Entropy_{modoList[1]}/entropy-{cList[1][i]:.3f}.txt', unpack=True)
    
    tlist_r, velocidade_r = Read_Velocidade(curvasList[0][i], modoList[0])
    tlist_a, velocidade_a = Read_Velocidade(curva, modoList[1])
    
    inter = np.argwhere(np.diff(np.sign(velocidade_a - velocidade_r))).flatten()
    
    tlist_inter = tlist[inter]
    va_inter = velocidade_a[inter]
    vr_inter = velocidade_r[inter]
    Sa_inter.append(S_a[inter])
    Sr_inter.append(S_r[inter])
    
    Sa_max.append(max(S_a))
    Sr_max.append(max(S_r))
    
    plt.scatter(tlist_a, velocidade_a, color='red', label='Heating')
    plt.scatter(tlist_r, velocidade_r, color='blue', label='Cooling')
    plt.scatter(tlist_inter, va_inter, color='red')
    plt.scatter(tlist_inter, vr_inter, color='blue')
    plt.xscale('log')
    plt.legend()
    plt.show()

plt.scatter(cList[0][:len(Sa_inter)], Sa_inter, color='red', label='Heating Intersection')
plt.scatter(cList[0][:len(Sr_inter)], Sr_inter, color='blue', label='Cooling Intersection')
plt.plot(cList[0][:len(Sa_max)], Sa_max, color='red', label='Heating Maximum')
plt.plot(cList[0][:len(Sr_max)], Sr_max, color='blue', label='Cooling Maximum')
plt.legend()
plt.show()