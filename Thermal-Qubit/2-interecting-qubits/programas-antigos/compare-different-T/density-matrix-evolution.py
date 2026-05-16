import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
import os


def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))


def Temperatura(w0, p):
    
    T = w0/np.log(p/(1-p))
    
    return T


def Entropia_Relativa_Populacoes(pi, pt):
    
    return pi*np.log(pi/pt) + (1-pi)*np.log((1-pi)/(1-pt))


####### MAIN #######

cmod = np.loadtxt(f'./DensityMatrices/cmod.txt', unpack=True)

cmod_extremes = [min(cmod), max(cmod)]

w0 = 1
Tf_qubit = 1.0
Sr_inicial = 0.05

tlist = np.arange(0, 20, 0.01)
timestep = np.arange(0, len(tlist), 1)

T = []
S = []
coherence_c = []

for c in cmod_extremes:
    
    print(c)
    
    Tt1 = []
    Tt2 = []
    Sc1 = []
    Sc2 = []
    coherence = []

    path = f'./DensityMatrices/c_{c}/'
    
    for t in timestep:
        
        rhot1 = np.loadtxt(path + f'rhof_q1_t{t}.txt', unpack=True, dtype='complex')
        rhot2 = np.loadtxt(path + f'rhof_q2_t{t}.txt', unpack=True, dtype='complex')
        rhoft = np.loadtxt(path + f'rhof_t{t}.txt', unpack=True, dtype='complex')
        
        p1 = rhot1[0][0].real
        p2 = rhot2[0][0].real

        T1 = Temperatura(w0, p1)
        Tt1.append(T1)
        
        T2 = Temperatura(w0, p2)
        Tt2.append(T2)
        
        Sc1.append( Entropia_Relativa_Populacoes(p1, pFunc(Tf_qubit, w0)) )
        Sc2.append( Entropia_Relativa_Populacoes(p2, pFunc(Tf_qubit, w0)) )
        
        fora_diagonal = rhoft[0][1] + rhoft[0][2] + rhoft[0][3] + rhoft[1][2] + rhoft[1][3] + rhoft[2][3]
        coherence.append( abs(fora_diagonal) )

    T.append(Tt1)
    T.append(Tt2)
    S.append(Sc1)
    S.append(Sc2)
    
    coherence_c.append(coherence)


plt.plot(tlist, T[0], label=f'Qubit 1')
plt.plot(tlist, T[1], label=f'Qubit 2')
plt.hlines(y=Tf_qubit, xmin=min(tlist), xmax=max(tlist), color='black', linestyle='--', label='Final Temperature')
plt.ylabel('Temperature')
plt.xlabel('Time')
plt.title(f'|c| = {cmod_extremes[0]:.3f}')
plt.legend()
plt.show()

plt.plot(tlist, T[2], label=f'Qubit 1')
plt.plot(tlist, T[3], label=f'Qubit 2')
plt.hlines(y=Tf_qubit, xmin=min(tlist), xmax=max(tlist), color='black', linestyle='--', label='Final Temperature')
plt.ylabel('Temperature')
plt.xlabel('Time')
plt.title(f'|c| = {cmod_extremes[1]:.3f}')
plt.legend()
plt.show()


plt.plot(tlist, S[0], label='Qubit 1')
plt.plot(tlist, S[1], label='Qubit 2')
plt.hlines(y=Sr_inicial, xmin=min(tlist), xmax=max(tlist), color='black', label='Initial Relative Entropy')
plt.ylabel('Relative Entropy')
plt.xlabel('Time')
plt.title(f'|c| = {cmod_extremes[0]:.3f}')
plt.legend()
plt.show()

plt.plot(tlist, S[2], label='Qubit 1')
plt.plot(tlist, S[3], label='Qubit 2')
plt.hlines(y=Sr_inicial, xmin=min(tlist), xmax=max(tlist), color='black', label='Initial Relative Entropy')
plt.ylabel('Relative Entropy')
plt.xlabel('Time')
plt.title(f'|c| = {cmod_extremes[1]:.3f}')
plt.legend()
plt.show()


plt.plot(tlist, coherence_c[0], label=f'|c| = {cmod_extremes[0]:.6f}')
plt.plot(tlist, coherence_c[1], label=f'|c| = {cmod_extremes[1]:.6f}')
plt.ylabel('Total Coherence')
plt.xlabel('Time')
plt.legend()
plt.show()
