import numpy as np
import matplotlib.pyplot as plt
import os


def Read_Completion(i, modo, completionList):

    path = f'./ThermalKinematics_{modo}/completion_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    t = np.where(0.21 == tlist)
    
    completionList.append(completion[t])
    
    return tlist, completionList
    
    
def Maximos(i, modo, maximos):

    curve_path = f'./FisherInformation_{modo}/curva_{i}.txt'
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    maximos.append(max(QFI))
    
    return maximos
    

def Coerencia(i, modo, cvalue):

    path = f'./FisherInformation_{modo}/c_curva_{i}.txt'
    
    c_all = np.loadtxt(path, unpack=True, dtype=complex, ndmin=1)
    
    cvalue.append(abs(c_all[0]))
    
    return cvalue
    

def Ordenamento(curvas, dados):

    dados, curvas_ordenadas = (list(t) for t in zip(*sorted(zip(dados, curvas))))
    
    return curvas_ordenadas

### MAIN ###


modoList = ['Resfriar', 'Aquecer']
maximos_modo = []
coerencia_modo = []

for modo in modoList:
    
    iList = []
    
    for arquivo in os.listdir(f'./ThermalKinematics_{modo}/'):
        
        if os.path.isfile(os.path.join(f'./ThermalKinematics_{modo}/', arquivo)) and 'completion_' in arquivo:
            
            iList.append(int(arquivo.replace('completion_', '').replace('.txt', '')))
    
    iList = sorted(iList)
    
    curva = []
    maximos = []
    cvalue = []
    completionList = []
    
    for i in iList:
    
        curva.append(i)
        
        tlist, completionList = Read_Completion(i, modo, completionList)
        
        maximos = Maximos(i, modo, maximos)

        cvalue = Coerencia(i, modo, cvalue)
        
    
    print(f'{modo}')
    
    print(f'Completion: {Ordenamento(curva, completionList)}')
    
    print(f'Maximos: {Ordenamento(curva, maximos)}')
    
    print(f'Coerência: {Ordenamento(curva, cvalue)}')
    
    print(' ')
    
    maximos_modo.append(sorted(maximos))
    coerencia_modo.append(sorted(cvalue))
    
plt.plot(coerencia_modo[0], maximos_modo[0], color='blue')
plt.scatter(coerencia_modo[0], maximos_modo[0], color='blue')
plt.xlabel('|c|')
plt.ylabel('QFI maximum')
plt.title('Cooling')
plt.tight_layout()
plt.show()

plt.plot(coerencia_modo[1], maximos_modo[1], color='red')
plt.scatter(coerencia_modo[1], maximos_modo[1], color='red')
plt.xlabel('|c|')
plt.ylabel('QFI maximum')
plt.title('Heating')
plt.tight_layout()
plt.show()


























