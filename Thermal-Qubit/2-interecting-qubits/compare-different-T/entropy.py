import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
import os
    

def RHO(timestep, c, path):
    
    rho_q1 = []
    rho_q2 = []
    
    for t in timestep:
        
        rhot1 = np.loadtxt(path + f'rhof_q1_t{t}.txt', unpack=True, dtype='complex')
        rhot2 = np.loadtxt(path + f'rhof_q2_t{t}.txt', unpack=True, dtype='complex')
        
        rho_q1.append(rhot1)
        rho_q2.append(rhot2)
	

    return rho_q1, rho_q2


def vonNeumann_Entropy(rho):
    
    S = []
    
    for rhot in rho:
        
        St = 0
        
        evals = Qobj(rhot).eigenenergies()

        for eval in evals:
            St = St - eval * np.log(eval)
        
        S.append(St)

    return S


def Derivada(xlist, ylist):
    
    ## https://www.youtube.com/watch?v=utRKIlOZbtw
    
    yprime = np.diff(ylist)/np.diff(xlist)
    xprime = []
    
    for i in range(len(yprime)):
        
        xtemp = (xlist[i+1] + xlist[i])/2
        xprime = np.append(xprime, xtemp)
    
    return xprime, yprime


def WriteFile(pasta, S, tlist):
    
    file = open(pasta, 'w')
    
    for k in range(len(tlist)):
        file.write(f'{tlist[k]} {S[k]}\n')
    
    file.close()
    

### MAIN ###

g = 0.8


tlist = np.arange(0.1, 20, 0.01)
timestep = np.arange(0, len(tlist), 1)

cmodlist = np.loadtxt(f'./g-{g}/DensityMatrices/cmod.txt', unpack=True)
cmodlist.sort()

os.mkdir(f'./g-{g}/Entropy')

Slist1 = []
Slist2 = []

for c in cmodlist:
    
    print(f'{c:.3f}')
    
    rho_q1, rho_q2 = RHO(timestep, c, f'./g-{g}/DensityMatrices/c_{c}/')
    
    S1 = vonNeumann_Entropy(rho_q1)
    S2 = vonNeumann_Entropy(rho_q2)
    
    equilibrio1 = np.where(S1==S1[-1])
    equilibrio2 = np.where(S2==S2[-1])
    
    print(f't={tlist[equilibrio1[0][0]]} ii={equilibrio1[0][0]} if={equilibrio1[0][-1]}')
    print(f't={tlist[equilibrio2[0][0]]} ii={equilibrio2[0][0]} if={equilibrio2[0][-1]}')
    
    Slist1.append(S1)
    Slist2.append(S2)
    
    WriteFile(f'./g-{g}/Entropy/entropy1-{c:.3f}.txt', S1, tlist)
    WriteFile(f'./g-{g}/Entropy/entropy2-{c:.3f}.txt', S2, tlist)
  
  

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

for i in range(len(cmodlist)):
    cor = next(colors)
    plt.plot(tlist, Slist1[i], color=cor, label=f'|c| = {cmodlist[i]:.3f}')    

plt.ylabel('S')
plt.xlabel('Time')
plt.title('Qubit 1')
plt.xscale('log')
plt.xlim(left=0.1)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()


colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

for i in range(len(cmodlist)):
    cor = next(colors)
    plt.plot(tlist, Slist2[i], color=cor, label=f'|c| = {cmodlist[i]:.3f}')    

plt.ylabel('S')
plt.xlabel('Time')
plt.title('Qubit 2')
plt.xscale('log')
plt.xlim(left=0.1)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()


colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

for s in range(len(Slist1)):
    
    dtlist, dS = Derivada(tlist, Slist1[s])
    
    WriteFile(f'./g-{g}/Entropy/d-entropy1-{cmodlist[s]:.3f}.txt', dS, dtlist)
    
    c = next(colors)
    plt.plot(dtlist, dS, color=c, label=f'|c| = {cmodlist[s]:.3f}')

plt.ylabel('dS/dt')
plt.xlabel('Time')
plt.title('Qubit 1')
plt.xscale('log')
plt.xlim(left=0.1)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()


colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

for s in range(len(Slist2)):
    
    dtlist, dS = Derivada(tlist, Slist2[s])
    
    WriteFile(f'./g-{g}/Entropy/d-entropy2-{cmodlist[s]:.3f}.txt', dS, dtlist)
    
    c = next(colors)
    plt.plot(dtlist, dS, color=c, label=f'|c| = {cmodlist[s]:.3f}')

plt.ylabel('dS/dt')
plt.xlabel('Time')
plt.title('Qubit 2')
plt.xscale('log')
plt.xlim(left=0.1)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()
