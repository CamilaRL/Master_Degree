import numpy as np
import matplotlib.pyplot as plt


def Integracao(ydata, xdata):

    dt = xdata[1] - xdata[0]
    
    L = np.cumsum(ydata) * dt

    return L


g = 0.8

cmod = np.loadtxt(f'./g-{g}/DensityMatrices/cmod.txt', unpack=True)

cmod_min = min(cmod)
cmod_max = max(cmod)


tQ1min, Q1min = np.loadtxt(f'./g-{g}/DensityMatrices/Q1_c{cmod_min}.txt', unpack=True)
tQ1max, Q1max = np.loadtxt(f'./g-{g}/DensityMatrices/Q1_c{cmod_max}.txt', unpack=True)

tQ2min, Q2min = np.loadtxt(f'./g-{g}/DensityMatrices/Q2_c{cmod_min}.txt', unpack=True)
tQ2max, Q2max = np.loadtxt(f'./g-{g}/DensityMatrices/Q2_c{cmod_max}.txt', unpack=True)

tSr1min, dSr1min = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q1_c{cmod_min}.txt', unpack=True)
tSr1max, dSr1max = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q1_c{cmod_max}.txt', unpack=True)

tSr2min, dSr2min = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q2_c{cmod_min}.txt', unpack=True)
tSr2max, dSr2max = np.loadtxt(f'./g-{g}/DensityMatrices/dSr_q2_c{cmod_max}.txt', unpack=True)


I_Q1min = Integracao(Q1min, tQ1min)
I_Q1max = Integracao(Q1max, tQ1max)

I_Q2min = Integracao(Q2min, tQ2min)
I_Q2max = Integracao(Q2max, tQ2max)

I_Sr1min = Integracao(dSr1min, tSr1min)
I_Sr1max = Integracao(dSr1max, tSr1max)

I_Sr2min = Integracao(dSr2min, tSr2min)
I_Sr2max = Integracao(dSr2max, tSr2max)


n_1min = []
n_1max = []
n_2min = []
n_2max = []

for t in range(len(tQ1min)):
    
    n_1min.append( I_Sr1min[t] / abs(I_Q1min[t]) )
    n_1max.append( I_Sr1max[t] / abs(I_Q1max[t]) )
    
    n_2min.append( I_Sr2min[t] / abs(I_Q2min[t]) )
    n_2max.append( I_Sr2max[t] / abs(I_Q2max[t]) )


fig = plt.figure(figsize=(10,5))

plt.subplot(121)

plt.plot(tQ1min, n_1min, color='red', label=f'|c| = {cmod_min:.3f}')
plt.plot(tQ1max, n_1max, color='orange', linestyle='--', label=f'|c| = {cmod_max:.3f}')

plt.title('Qubit 1')
plt.ylabel('Efficiency')
plt.xlabel('Time')
plt.legend()

plt.subplot(122)

plt.plot(tQ2min, n_2min, color='blue', label=f'|c| = {cmod_min:.3f}')
plt.plot(tQ2max, n_2max, color='skyblue', linestyle='--', label=f'|c| = {cmod_max:.3f}')

plt.title('Qubit 2')
plt.ylabel('Efficiency')
plt.xlabel('Time')
plt.legend()


plt.tight_layout()
plt.show()

