import numpy as np
import matplotlib.pyplot as plt
import os


def Read_Completion(i, modo, completion1, completion2):

    path = f'./ThermalKinematics_{modo}/completion_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    volta1 = np.where(0.71 == tlist)
    volta2 = np.where(2.49 == tlist)
    
    completion1.append(completion[volta1])
    completion2.append(completion[volta2])
    
    return tlist, completion1, completion2
    
    
def Maximos(i, modo, maximos):

    curve_path = f'./FisherInformation_{modo}/curva_{i}.txt'
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    for t in tlist:
    
        if np.sin(2*t) == 0 and np.cos(2*t) == -1:
            
            i = np.where(tlist == t)
            break
    
    maximos.append(QFI[i])
    
    #maximos.append(max(QFI))
    
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
    completion1 = []
    completion2 = []
    
    for i in iList:
    
        curva.append(i)
        
        tlist, completion1, completion2 = Read_Completion(i, modo, completion1, completion2)
        
        maximos = Maximos(i, modo, maximos)

        cvalue = Coerencia(i, modo, cvalue)
        
    
    print(f'{modo}')
    
    print(f'Completion volta 1: {Ordenamento(curva, completion1)}')
   
    print(f'Completion volta 2: {Ordenamento(curva, completion2)}')
    
    print(f'Maximos: {Ordenamento(curva, maximos)}')
    
    print(f'Coerência: {Ordenamento(curva, cvalue)}')
    
    print(' ')
    
    maximos_modo.append(sorted(maximos))
    coerencia_modo.append(sorted(cvalue))
    
plt.plot(coerencia_modo[0], maximos_modo[0], color='blue', label=modoList[0])
plt.scatter(coerencia_modo[0], maximos_modo[0], color='blue')
plt.plot(coerencia_modo[1], maximos_modo[1], color='red', label=modoList[1])
plt.scatter(coerencia_modo[1], maximos_modo[1], color='red')

plt.xlabel('Amplitude de Coerência')
plt.ylabel('Maximos de QFI')
plt.legend()
plt.tight_layout()
plt.show()

























