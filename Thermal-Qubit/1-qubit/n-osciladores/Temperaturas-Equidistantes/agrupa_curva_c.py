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


dTList = [1, 5 ,9]
modoList = ['Resfriar', 'Aquecer']

cmodList = []
curvasList = []


for modo in modoList:

    matrix_agrupamento = []

    for dT in dTList:

        iList = []
        cvalue = []
        
        for arquivo in os.listdir(f'./dT-{dT}/FisherInformation_{modo}/'):

            if os.path.isfile(os.path.join(f'./dT-{dT}/FisherInformation_{modo}/', arquivo)) and 'c_curva_' in arquivo:
                
                curva_i = int(arquivo.replace('c_curva_', '').replace('.txt', ''))
                
                if curva_i > 0:
                
                    iList.append(curva_i)

        for i in iList:
            
            cvalue.append(Read_Mod_Cvalue(i, dT, modo))
        
        cmod, iList = Ordenamento(cvalue, iList)
        
        if len(matrix_agrupamento) < 1:
        
            matrix_agrupamento.append(cmod)
            matrix_agrupamento.append(iList)
            
        else:
        
            nova_coluna = [0 for i in range(len(matrix_agrupamento[0]))]

            for j, c in enumerate(cmod):
                
                if c in matrix_agrupamento[0]:
                    
                    k = matrix_agrupamento[0].index(c)
                    
                    nova_coluna[k] = iList[j]
                    
                else:
                
                    matrix_agrupamento[0].append(c)
                    for m in range(1, len(matrix_agrupamento), 1):
                    
                        matrix_agrupamento[m].append(0)
                        
                    nova_coluna.append(iList[j])
            
            matrix_agrupamento.append(nova_coluna)
    
    f_a = open(f'./Agrupamento/agrupamento_{modo}.txt', 'w')
    
    for i in range(len(matrix_agrupamento[0])):
        f_a.write(f'{matrix_agrupamento[0][i]} {matrix_agrupamento[1][i]} {matrix_agrupamento[2][i]} {matrix_agrupamento[3][i]}\n')
    
    f_a.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    
