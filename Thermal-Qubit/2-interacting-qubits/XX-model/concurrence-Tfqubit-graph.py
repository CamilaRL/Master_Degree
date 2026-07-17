import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
import os


def vonNeumann_Entropy(rho):
    
    St = 0
        
    evals = rho.eigenenergies()
        
    for ev in evals:
        if ev > 0:
            St = St - ev * np.log(ev)
        
    return St
    

w0 = 1
Jlist = np.arange(0, 10, 0.01)
beta_R = 1

Clist = []
Elist = []
MIlist = []
betaList = []

for J in Jlist:

    ## concurrence
    
    Z = 2*(np.cosh(beta_R*w0) + np.cosh(beta_R*J))

    lambdas = (np.exp(beta_R*J) - 2 - np.exp(-beta_R*J))/Z
    
    C = max([0, lambdas])
    
    Clist.append(C) 
    
    
    ## EoF
    
    x = 0.5*(1 + np.sqrt(1 - C**2))

    if x == 0 or x == 1:
        E = 0
            
    else:
        E = -x * math.log2(x) - (1-x) * math.log2(1-x)
        
    Elist.append(E)
    
    
    ## Mutual Information
    
    rhoq = Qobj([[(1/Z)*(np.exp(beta_R*w0)+np.cosh(beta_R*J)), 0], [0, (1/Z)*(np.exp(-beta_R*w0)+np.cosh(beta_R*J))]])
    rhoG = Qobj([[(1/Z)*np.exp(beta_R*w0), 0, 0, 0], [0, (1/Z)*np.cosh(beta_R*J), -(1/Z)*np.sinh(beta_R*J), 0], [0, -(1/Z)*np.sinh(beta_R*J), (1/Z)*np.cosh(beta_R*J), 0], [0, 0, 0, (1/Z)*np.exp(-beta_R*w0)]])
    
    S1 = vonNeumann_Entropy(rhoq)
    S2 = vonNeumann_Entropy(rhoq)
    S12 = vonNeumann_Entropy(rhoG)
    
    MIlist.append(S1 + S2 - S12)
    
    
    ## Temperatura final dos qubits
    
    beta = (1/w0)*np.log((np.exp(beta_R*w0) + np.cosh(beta_R*J))/(np.exp(-beta_R*w0) + np.cosh(beta_R*J)))
    
    betaList.append(beta)
    
    
### PLOTS

fig = plt.figure(figsize=(12,5))

plt.subplot(131)
plt.plot(Jlist, Clist, color='black', linewidth=2)
plt.xlabel(r'$J$', fontsize=12)
plt.ylabel('Concurrence', fontsize=12)

plt.subplot(132)
plt.plot(Jlist, Elist, color='black', linewidth=2)
plt.xlabel(r'$J$', fontsize=12)
plt.ylabel('Entanglement of Formation', fontsize=12)

plt.subplot(133)
plt.plot(Jlist, MIlist, color='black', linewidth=2)
plt.xlabel(r'$J$', fontsize=12)
plt.ylabel('Mutual Information', fontsize=12)

plt.tight_layout()
plt.show()


plt.plot(Jlist, betaList, color='black', linewidth=2)
plt.xlabel(r'$J$', fontsize=12)
plt.ylabel(r'$\beta_{f}^{q}$', fontsize=12)
plt.show()