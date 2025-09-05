import numpy as np
import matplotlib.pyplot as plt
import os


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
modoIngles = ['Cooling', 'Heating']
cores = ['blue', 'red']


for modo in modoList:
    
    m = modoList.index(modo)
    
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


    t_equilibrio = []
    cmod = []
    
    cmap = plt.get_cmap('rainbow')
    colors = iter(cmap(np.linspace(0.01, 1, len(iList))))


    for i in range(len(iList)):

        tlist, completion = np.loadtxt(f'./ThermalKinematics_{modo}/completion_{iList[i]}.txt', unpack=True)        
        
        equilibrio = np.where(completion==completion[-1])
        
        t_equilibrio.append(tlist[equilibrio[0][0]])
        
        cmod.append(abs(cvalue[i]))
        
        

    plt.scatter(cmod, t_equilibrio, color=cores[m])
    plt.title(f'{modoIngles[m]}')
    plt.ylabel('Equilibration Time')
    plt.xlabel('|c|')
    plt.tight_layout()
    plt.show()
