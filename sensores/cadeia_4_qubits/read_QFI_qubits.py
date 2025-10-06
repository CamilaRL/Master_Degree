import numpy as np
import matplotlib.pyplot as plt

gList = [0.01, 0.5, 0.8, 1.0, 1.5]
JList = [1.0, 2.0]

fig = plt.figure(figsize=(10,8))

for J in JList:

    tSEmax = 10*J
    
    QFI_max1 = []
    QFI_max2 = []
    QFI_max3 = []
    QFI_max4 = []
    
    for g in gList:

        temp1, QFI_q1 = np.loadtxt(f'./Results/QFI_1_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp2, QFI_q2 = np.loadtxt(f'./Results/QFI_2_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp3, QFI_q3 = np.loadtxt(f'./Results/QFI_3_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
        temp4, QFI_q4 = np.loadtxt(f'./Results/QFI_4_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)

        QFI_max1.append(max(QFI_q1))
        QFI_max2.append(max(QFI_q2))
        QFI_max3.append(max(QFI_q3))
        QFI_max4.append(max(QFI_q4))
    
    plt.subplot(221)
    plt.scatter(gList, QFI_max1, s=10, label=f'J = {J}')
    plt.plot(gList, QFI_max1)
    
    plt.subplot(222)
    plt.scatter(gList, QFI_max2, s=10, label=f'J = {J}')
    plt.plot(gList, QFI_max2)
    
    plt.subplot(223)
    plt.scatter(gList, QFI_max3, s=10, label=f'J = {J}')
    plt.plot(gList, QFI_max3)
    
    plt.subplot(224)
    plt.scatter(gList, QFI_max4, s=10, label=f'J = {J}')
    plt.plot(gList, QFI_max4)
    
    
plt.subplot(221)
plt.title('Qubit 1')
plt.ylabel('Maximum QFI')
plt.legend()

plt.subplot(222)
plt.title('Qubit 2')
plt.legend()

plt.subplot(223)
plt.title('Qubit 3')
plt.ylabel('Maximum QFI')
plt.xlabel('g')
plt.legend()

plt.subplot(224)
plt.title('Qubit 4')
plt.xlabel('g')
plt.legend()

plt.tight_layout()
plt.show()


