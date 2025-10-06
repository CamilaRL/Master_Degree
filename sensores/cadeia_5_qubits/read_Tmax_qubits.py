import numpy as np
import matplotlib.pyplot as plt

gList = [0.01, 0.5, 1.0, 1.5]
JList = [1.0, 2.0]

fig = plt.figure(figsize=(10,8))

for J in JList:

    tSEmax = 10*J
    
    T_max1 = []
    T_max2 = []
    T_max3 = []
    T_max4 = []
    T_max5 = []
    
    for g in gList:

        temp1, QFI_q1 = np.loadtxt(f'./Results/QFI_1_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp2, QFI_q2 = np.loadtxt(f'./Results/QFI_2_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp3, QFI_q3 = np.loadtxt(f'./Results/QFI_3_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp4, QFI_q4 = np.loadtxt(f'./Results/QFI_4_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp5, QFI_q5 = np.loadtxt(f'./Results/QFI_5_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        
        T_max1.append(temp1[np.argmax(QFI_q1)])
        T_max2.append(temp2[np.argmax(QFI_q2)])
        T_max3.append(temp3[np.argmax(QFI_q3)])
        T_max4.append(temp4[np.argmax(QFI_q4)])
        T_max5.append(temp5[np.argmax(QFI_q5)])
    
    
    plt.subplot(231)
    plt.scatter(gList, T_max1, s=10, label=f'J = {J}')
    plt.plot(gList, T_max1)
    
    plt.subplot(232)
    plt.scatter(gList, T_max2, s=10, label=f'J = {J}')
    plt.plot(gList, T_max2)
    
    plt.subplot(233)
    plt.scatter(gList, T_max3, s=10, label=f'J = {J}')
    plt.plot(gList, T_max3)
    
    plt.subplot(234)
    plt.scatter(gList, T_max4, s=10, label=f'J = {J}')
    plt.plot(gList, T_max4)
    
    plt.subplot(235)
    plt.scatter(gList, T_max5, s=10, label=f'J = {J}')
    plt.plot(gList, T_max5)
    
    
plt.subplot(231)
plt.title('Qubit 1')
plt.ylabel('Temperature of Maximum QFI')
plt.legend()

plt.subplot(232)
plt.title('Qubit 2')
plt.legend()

plt.subplot(233)
plt.title('Qubit 3')
plt.legend()

plt.subplot(234)
plt.ylabel('Temperature of Maximum QFI')
plt.title('Qubit 4')
plt.xlabel('g')
plt.legend()

plt.subplot(235)
plt.title('Qubit 5')
plt.xlabel('g')
plt.legend()

plt.tight_layout()
plt.show()

