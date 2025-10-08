import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
import os

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) + 1)
    

def RHO(tlist, c, p, gamma, w, nbar):
	
	rho = []
	
	for t in tlist:

		rx = c.real * np.exp(-2*gamma*(nbar + 0.5)*t)
		ry = -c.imag * np.exp(-2*gamma*(nbar + 0.5)*t)
		rz = (1/(2*nbar + 1)) - 2*((nbar + 1)/(2*nbar + 1) - p) * np.exp(-2*gamma*(2*nbar + 1)*t)

		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		rho.append(r_list)
	

	return rho


def vonNeumann_Entropy(rho):
    
    S = []
    
    for t in range(len(rho)):
        
        St = 0
        
        rhot = 0.5 * Qobj( [[1 + rho[t][2], complex(rho[t][0], -rho[t][1])], [complex(rho[t][0], rho[t][1]), 1 - rho[t][2]]] )
        
        evals = rhot.eigenenergies()

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

Sr = 0.1

modo = 'Cooling'

Tbanho = 1.4426950408889627
w = 2
Tqubit = 4.700782656252254
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))

gamma = 3

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)

coerencia_max = np.sqrt(max(0.0, p*(1-p)))

curvas = np.loadtxt(f'./FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)
cmodlist = np.loadtxt(f'./FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)

os.mkdir(f'./Entropy_{modo}_{Sr}')

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

Slist = []

for i in range(len(cmodlist)):
    
    print(cmodlist[i])
    
    
    rho = RHO(tlist, 0, p, gamma, w, nbar)
    
    S = vonNeumann_Entropy(rho)
    
    equilibrio = np.where(S==S[-1])
    
    print(f'c = 0 : t={tlist[equilibrio[0][0]]} ii={equilibrio[0][0]} if={equilibrio[0][-1]}')
    
    Slist.append(S)
    
    c = next(colors)
    
    WriteFile(f'./Entropy_{modo}_{Sr}/entropy-{cmodlist[i]:.3f}.txt', S, tlist)
    
    plt.plot(tlist, S, color=c, label=f'|c| = {cmodlist[i]:.3f}')

plt.ylabel('S')
plt.xlabel('Time')
plt.title(modo)
plt.xscale('log')
plt.xlim(left=0.01)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()

colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

for s in range(len(Slist)):
    
    dtlist, dS = Derivada(tlist, Slist[s])
    
    c = next(colors)
    plt.plot(dtlist, dS, color=c, label=f'|c| = {cmodlist[s]:.3f}')

plt.ylabel('dS/dt')
plt.xlabel('Time')
plt.title(modo)
plt.xscale('log')
plt.xlim(left=0.01)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()
