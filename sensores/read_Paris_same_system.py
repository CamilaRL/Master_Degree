import numpy as np
import matplotlib.pyplot as plt



g = 0.45
J = 1.0
tSEmax = 10*J

temp, QFI_all = np.loadtxt(f'./Results/QFI_all_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)

temp1, QFI_q1 = np.loadtxt(f'./Results/QFI_q1_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp2, QFI_q2 = np.loadtxt(f'./Results/QFI_q2_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp3, QFI_q3 = np.loadtxt(f'./Results/QFI_q3_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp4, QFI_q4 = np.loadtxt(f'./Results/QFI_q4_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)

'''QFI_all = []
for i in range(len(QFI_q1)):
    QFI_all.append( QFI_q1[i] + QFI_q2[i] + QFI_q3[i] + QFI_q4[i] )
'''

plt.plot(temp1, QFI_all, color='black', label='Total QFI')

plt.plot(temp1, QFI_q1, color='blue', linestyle='-', label=f'Qubit 1 - Max = {max(QFI_q1):.3f}')
plt.plot(temp2, QFI_q2, color='green', linestyle='-', label=f'Qubit 2 - Max = {max(QFI_q2):.3f}')
plt.plot(temp3, QFI_q3, color='orange', linestyle='--', label=f'Qubit 3 - Max = {max(QFI_q3):.3f}')
plt.plot(temp4, QFI_q4, color='red', linestyle='--', label=f'Qubit 4 - Max = {max(QFI_q4):.3f}')

plt.ylabel('QFI')
plt.xlabel('Temperature')
plt.legend()
plt.title(f'Transverse Field: {g:.2f} | Thermalization Time: {tSEmax:.1f}')
plt.show()
