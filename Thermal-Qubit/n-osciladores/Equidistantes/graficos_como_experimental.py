import numpy as np
import matplotlib.pyplot as plt


def Read_Completion(i, dT, modo):

    path = f'./dT-{dT}/ThermalKinematics_{modo}/completion_{i}.txt'

    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(i, dT, modo):

    path = f'./dT-{dT}/ThermalKinematics_{modo}/velocity_{i}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion



### MAIN ###

dTList = [1, 5, 9]
linhas = [':', '--', '-']
modoList = ['Resfriar', 'Aquecer']
modoIngles = ['Cooling', 'Heating']
cores = ['blue', 'red']

dT = []

cList_r, dT1, dT5, dT9 = np.loadtxt(f'./Agrupamento/agrupamento_{modoList[0]}.txt', unpack=True)        
dT.append([dT1, dT5, dT9])

cList_a, dT1, dT5, dT9 = np.loadtxt(f'./Agrupamento/agrupamento_{modoList[1]}.txt', unpack=True)
dT.append([dT1, dT5, dT9])

cList = []
        
for j, c in enumerate(cList_r):

    if c in cList_a:
        i = np.where(c == cList_a)[0][0]
        cList.append([c, j, i])

for citem in cList:

    c = citem[0]
    ir = citem[1]
    ia = citem[2]

    fig = plt.figure(figsize=(10,5))
    
    for d in range(3):

        curva_r = int(dT[0][d][ir])
        curva_a = int(dT[1][d][ia])
        
        
        if curva_r != 0 and curva_a != 0:
            
            legenda = [f'{modoIngles[0]} - ' + r'$\Delta$T = ' + f'{dTList[d]}', f'{modoIngles[1]} - ' + r'$\Delta$T = ' + f'{dTList[d]}']
                
            tlist_r, completion_r = Read_Completion(curva_r, dTList[d], modoList[0])
            tlist_a, completion_a = Read_Completion(curva_a, dTList[d], modoList[1])
                
            plt.subplot(121)
            plt.plot(tlist_r, completion_r, linestyle=linhas[d], color=cores[0], label=legenda[0])
            plt.plot(tlist_a, completion_a, linestyle=linhas[d], color=cores[1], label=legenda[1])
                
            
            tlist_r, velocidade_r = Read_Velocidade(curva_r, dTList[d], modoList[0])
            tlist_a, velocidade_a = Read_Velocidade(curva_a, dTList[d], modoList[1])
                
            plt.subplot(122)
            plt.plot(tlist_r, velocidade_r, linestyle=linhas[d], color=cores[0], label=legenda[0])
            plt.plot(tlist_a, velocidade_a, linestyle=linhas[d], color=cores[1], label=legenda[1])
            


    plt.subplot(121)
    plt.xlabel('Time')
    plt.ylabel('Degree of completion')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.ylim(bottom=0.5)


    plt.subplot(122)
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.legend(loc='center right')
    plt.xscale('log')
    plt.xlim(left=0.01)
    plt.ylim(top=40)
    


    plt.suptitle(f'|c| = {c:.3f}')
    plt.tight_layout()
    plt.show()
