import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
import os


def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))
    
    
def Entropia_Relativa_Populacoes(pi, pt):
    
    return pi*np.log(pi/pt) + (1-pi)*np.log((1-pi)/(1-pt))
    
    
def Derivada(xlist, ylist):
    
    ## https://www.youtube.com/watch?v=utRKIlOZbtw
    
    yprime = np.diff(ylist)/np.diff(xlist)
    xprime = []
    
    for i in range(len(yprime)):
        
        xtemp = (xlist[i+1] + xlist[i])/2
        xprime = np.append(xprime, xtemp)
    
    return xprime, yprime

    
  

######### MAIN #########

g = 0
w0 = 1
Tf_qubit = 1.0 #1.1542306681364487
Sr_inicial = 0.05


cmod = np.loadtxt(f'./g-{g}/DensityMatrices/cmod.txt', unpack=True)

cmod_extremes = [min(cmod), max(cmod)]

tlist = np.arange(0, 20, 0.01)
timestep = np.arange(0, len(tlist), 1)

for c in cmod_extremes:
    
    path = f'./g-{g}/DensityMatrices/'
    
    Sc1 = []
    Sc2 = []
    
    for t in timestep:
        
        rhot1 = np.loadtxt(path + f'c_{c}/rhof_q1_t{t}.txt', unpack=True, dtype='complex')
        rhot2 = np.loadtxt(path + f'c_{c}/rhof_q2_t{t}.txt', unpack=True, dtype='complex')
        
        p1 = rhot1[0][0].real
        p2 = rhot2[0][0].real
        
        Sc1.append( Entropia_Relativa_Populacoes(p1, pFunc(Tf_qubit, w0)) )
        Sc2.append( Entropia_Relativa_Populacoes(p2, pFunc(Tf_qubit, w0)) )
        
    dtlist1, dSc1 = Derivada(tlist, Sc1)
    dtlist2, dSc2 = Derivada(tlist, Sc2)
    
    
    for i in range(len(dSc1)):
        
        dSc1[i] = - dSc1[i]
        dSc2[i] = - dSc2[i]
    
    
    plt.plot(tlist, Sc1, label='Qubit 1')
    plt.plot(tlist, Sc2, label='Qubit 2')
    plt.ylabel('Relative Entropy')
    plt.xlabel('Time')
    plt.title(f'|c| = {c:.3f}')
    plt.legend()
    plt.show()
    
    plt.plot(dtlist1, dSc1, label='Qubit 1')
    plt.plot(dtlist2, dSc2, label='Qubit 2')
    plt.ylabel('Local Entropy Production Rate')
    plt.xlabel('Time')
    plt.title(f'|c| = {c:.3f}')
    plt.legend()
    plt.show()

    
    fSc1 = open(path + f'Sr_q1_c{c}.txt', 'w')
    fSc2 = open(path + f'Sr_q2_c{c}.txt', 'w')
    fdSc1 = open(path + f'dSr_q1_c{c}.txt', 'w')
    fdSc2 = open(path + f'dSr_q2_c{c}.txt', 'w')
    
    for t in range(len(tlist)):
        
        fSc1.write(f'{tlist[t]} {Sc1[t]}\n')
        fSc2.write(f'{tlist[t]} {Sc2[t]}\n')
        
    for t in range(len(dtlist1)):
        
        fdSc1.write(f'{dtlist1[t]} {dSc1[t]}\n')
        fdSc2.write(f'{dtlist2[t]} {dSc2[t]}\n')
