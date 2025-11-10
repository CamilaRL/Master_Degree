import numpy as np
import matplotlib.pyplot as plt


tot = 7

curva = []
cvalue = []
maximos = []
c_allList = []

fig = plt.figure(figsize=(10,5))

for i in range(tot):
    
    path = f'./FisherInformation_Resfriar/c_curva_{i}.txt'
    curve_path = f'./FisherInformation_Resfriar/curva_{i}.txt'
    
    curva.append(i)
    
    c_all = np.loadtxt(path, unpack=True, dtype=complex, ndmin=1)
    
    c_allList.append(c_all)
    
    tlist, QFI = np.loadtxt(curve_path, unpack=True)
    
    maximos.append(max(QFI))
    
    xlist = []
    clist = []
    
    for j in range(len(c_all)):
        
        xlist.append(i)
        clist.append(abs(c_all[j]))

    cvalue.append(clist[0])
    
    plt.subplot(122)
    plt.scatter(xlist, clist)




plt.subplot(121)
plt.plot(curva, maximos)
plt.scatter(curva, maximos)
plt.xlabel('Curvas')
plt.ylabel('Maximos')

plt.subplot(122)
plt.xlabel('Curvas')
plt.ylabel('Módulo da amplitude de coerência (c)')

plt.tight_layout()
plt.show()



