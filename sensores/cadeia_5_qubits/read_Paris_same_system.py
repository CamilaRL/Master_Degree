import numpy as np
import matplotlib.pyplot as plt



g = 1.5
J = 1.0
tSEmax = 10*J

temp, QFI_all = np.loadtxt(f'./Results/QFI_all_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)

temp1, QFI_q1 = np.loadtxt(f'./Results/QFI_1_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp2, QFI_q2 = np.loadtxt(f'./Results/QFI_2_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp3, QFI_q3 = np.loadtxt(f'./Results/QFI_3_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp4, QFI_q4 = np.loadtxt(f'./Results/QFI_4_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp5, QFI_q5 = np.loadtxt(f'./Results/QFI_5_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)

plt.plot(temp1, QFI_all, color='black', label='Total QFI')

plt.plot(temp1, QFI_q1, color='blue', linestyle='-', label=f'Qubit 1 - QFI_max = {max(QFI_q1):.3f}')
plt.plot(temp2, QFI_q2, color='green', linestyle='-', label=f'Qubit 2 - QFI_max = {max(QFI_q2):.3f}')
plt.plot(temp3, QFI_q3, color='orange', linestyle='--', label=f'Qubit 3 - QFI_max = {max(QFI_q3):.3f}')
plt.plot(temp4, QFI_q4, color='red', linestyle=':', label=f'Qubit 4 - QFI_max = {max(QFI_q4):.3f}')
plt.plot(temp5, QFI_q5, color='purple', linestyle='--', label=f'Qubit 5 - QFI_max = {max(QFI_q5):.3f}')

plt.ylabel('QFI')
plt.xlabel('Temperature')
plt.legend()
plt.title(f'Transverse Field: {g:.2f} | Thermalization Time: {tSEmax:.1f}')
plt.show()
