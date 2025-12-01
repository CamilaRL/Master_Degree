import numpy as np
import matplotlib.pyplot as plt



g = 0.8

cmod = np.loadtxt(f'./g-{g}/DensityMatrices/cmod.txt', unpack=True)

cmod_min = min(cmod)
cmod_max = max(cmod)


##### relative entropy

tlist1_c0, Sr1_c0 = np.loadtxt(f'./g-{g}/DensityMatrices/Sr_q1_c{cmod_min}.txt', unpack=True)
tlist2_c0, Sr2_c0 = np.loadtxt(f'./g-{g}/DensityMatrices/Sr_q2_c{cmod_min}.txt', unpack=True)

tlist1_cmax, Sr1_cmax = np.loadtxt(f'./g-{g}/DensityMatrices/Sr_q1_c{cmod_max}.txt', unpack=True)
tlist2_cmax, Sr2_cmax = np.loadtxt(f'./g-{g}/DensityMatrices/Sr_q2_c{cmod_max}.txt', unpack=True)
    
    
plt.plot(tlist1_c0, Sr1_c0, color='red', linestyle='-', label=f'|c| = {cmod_min:.3f} - Qubit 1')
plt.plot(tlist2_c0, Sr2_c0, color='blue', linestyle='-', label=f'|c| = {cmod_min:.3f} - Qubit 2')

plt.plot(tlist1_cmax, Sr1_cmax, color='orange', linestyle='--', label=f'|c| = {cmod_max:.3f} - Qubit 1')
plt.plot(tlist2_cmax, Sr2_cmax, color='skyblue', linestyle='--', label=f'|c| = {cmod_max:.3f} - Qubit 2')

plt.ylabel('Relative Entropy', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title(f'g = {g}')
plt.legend(fontsize=12)
plt.show()


##### entropy production

dtlist1_c0, dSr1_c0 = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q1_c{cmod_min}.txt', unpack=True)
dtlist2_c0, dSr2_c0 = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q2_c{cmod_min}.txt', unpack=True)

dtlist1_cmax, dSr1_cmax = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q1_c{cmod_max}.txt', unpack=True)
dtlist2_cmax, dSr2_cmax = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q2_c{cmod_max}.txt', unpack=True)


    
plt.plot(dtlist1_c0, dSr1_c0, color='red', linestyle='-', linewidth=2, label=f'|c| = {cmod_min:.3f} - Qubit 1')
plt.plot(dtlist2_c0, dSr2_c0, color='blue', linestyle='-', linewidth=2, label=f'|c| = {cmod_min:.3f} - Qubit 2')

plt.plot(dtlist1_cmax, dSr1_cmax, color='orange', linestyle='-', linewidth=2, label=f'|c| = {cmod_max:.3f} - Qubit 1')
plt.plot(dtlist2_cmax, dSr2_cmax, color='skyblue', linestyle='-', linewidth=2, label=f'|c| = {cmod_max:.3f} - Qubit 2')

plt.ylabel('Local Entropy Production Rate', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title(f'g = {g}')
plt.legend(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.show()