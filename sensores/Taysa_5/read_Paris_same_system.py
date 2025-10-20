import numpy as np
import matplotlib.pyplot as plt



g = 0.45
J = 1.0
tSEmax = 10*J

temp, QFI_all = np.loadtxt(f'./Results/QFI_all_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)

temp1, QFI_q1 = np.loadtxt(f'./Results/QFI_q1_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp2, QFI_q2 = np.loadtxt(f'./Results/QFI_q2_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp3, QFI_q3 = np.loadtxt(f'./Results/QFI_q3_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp4, QFI_q4 = np.loadtxt(f'./Results/QFI_q4_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp5, QFI_q5 = np.loadtxt(f'./Results/QFI_q5_g{g:.2f}_ttherm{tSEmax:.3f}.txt', unpack=True)

plt.plot(temp, QFI_all, color='black', label='Total QFI')

plt.plot(temp1, QFI_q1, color='blue', linestyle='-', label=f'Qubit 1 - TMax = {temp1[np.argmax(QFI_q1)]:.3f}')
plt.plot(temp2, QFI_q2, color='green', linestyle='-', label=f'Qubit 2 - TMax = {temp2[np.argmax(QFI_q1)]:.3f}')
plt.plot(temp3, QFI_q3, color='orange', linestyle='--', label=f'Qubit 3 - TMax = {temp3[np.argmax(QFI_q1)]:.3f}')
plt.plot(temp4, QFI_q4, color='red', linestyle='--', label=f'Qubit 4 - TMax = {temp4[np.argmax(QFI_q1)]:.3f}')
plt.plot(temp5, QFI_q5, color='purple', linestyle='--', label=f'Qubit 5 - TMax = {temp5[np.argmax(QFI_q1)]:.3f}')

plt.ylabel('QFI')
plt.xlabel('Temperature')
plt.legend()
plt.title(f'Transverse Field: {g:.2f} | Thermalization Time: {tSEmax:.1f}')
plt.show()
