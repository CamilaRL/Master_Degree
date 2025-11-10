import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

def func_Completion_dif(x, a, b):
    return a*x + b

def func_Velocidade_dif(x, a, b, c):
    return a*np.exp(b*x) + c


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


completion_dif = []
velocidade_dif = []

for i, curva in enumerate(curvasList[1]):
    
    tlist_r, completion_r = Read_Completion(curvasList[0][i], modoList[0])
    tlist_a, completion_a = Read_Completion(curva, modoList[1])
    
    completion_dif.append(completion_a[1] - completion_r[1])
    
    tlist_r, velocidade_r = Read_Velocidade(curvasList[0][i], modoList[0])
    tlist_a, velocidade_a = Read_Velocidade(curva, modoList[1])
    
    velocidade_dif.append(velocidade_a[1] - velocidade_r[1])
    

popt_completion, pcov_completion = curve_fit(func_Completion_dif, cList[1][:-2], completion_dif[:-2])
popt_velocidade, pcov_velocidade = curve_fit(func_Velocidade_dif, cList[1][:-2], velocidade_dif[:-2])

completion_fit = []
velocidade_fit = []
for c in cList[1][:-2]:
    completion_fit.append(func_Completion_dif(c, *popt_completion))
    velocidade_fit.append(func_Velocidade_dif(c, *popt_velocidade))
    

plt.scatter(cList[1][:-2], completion_dif[:-2])
plt.plot(cList[1][:-2], completion_fit)
plt.title('Difference of Initial Degree of Completion')
plt.ylabel('Degree of Completion Difference')
plt.xlabel('|c|')
plt.show()

plt.scatter(cList[1][:-2], velocidade_dif[:-2])
plt.plot(cList[1][:-2], velocidade_fit)
plt.title('Difference of Initial Velocity')
plt.ylabel('Velocity Difference')
plt.xlabel('|c|')
plt.show()
    