import numpy as np
import matplotlib.pyplot as plt



g = 0.45
tSEmax = (np.pi/2)+100

temp, QFI_all = np.loadtxt(f'./Results/k4QFI_all_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)

temp1, QFI_q1 = np.loadtxt(f'./Results/k4QFI_q1_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp2, QFI_q2 = np.loadtxt(f'./Results/k4QFI_q2_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp3, QFI_q3 = np.loadtxt(f'./Results/k4QFI_q3_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)
temp4, QFI_q4 = np.loadtxt(f'./Results/k4QFI_q4_g{g:.1f}_ttherm{tSEmax:.3f}.txt', unpack=True)


temp_max_all = temp[np.argmax(QFI_all)]
temp_max_q1 = temp[np.argmax(QFI_q1)]
temp_max_q2 = temp[np.argmax(QFI_q2)]
temp_max_q3 = temp[np.argmax(QFI_q3)]
temp_max_q4 = temp[np.argmax(QFI_q4)]


plt.plot(temp, QFI_all, color='black', label=f'Total QFI - T_max = {temp_max_all:.3f}')

plt.plot(temp1, QFI_q1, color='blue', linestyle='-', label=f'Qubit 1 - T_max = {temp_max_q1:.3f}')
plt.plot(temp2, QFI_q2, color='green', linestyle='-', label=f'Qubit 2 - T_max = {temp_max_q2:.3f}')
plt.plot(temp3, QFI_q3, color='orange', linestyle='--', label=f'Qubit 3 - T_max = {temp_max_q3:.3f}')
plt.plot(temp4, QFI_q4, color='red', linestyle='--', label=f'Qubit 4 - T_max = {temp_max_q4:.3f}')

plt.ylabel('QFI')
plt.xlabel('Temperature')
plt.legend()
plt.title(f'Transverse Field: {g:.2f} | Thermalization Time: {tSEmax:.3f}')
plt.show()
