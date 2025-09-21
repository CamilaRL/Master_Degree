import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) - 1)
    

def RHO(tlist, c, p, gamma, w, nbar):
	
	rho = []
	rho_derivada = []
	
	for t in tlist:

		rx = c.real * np.exp(-2*gamma*(nbar + 0.5)*t)
		ry = -c.imag * np.exp(-2*gamma*(nbar + 0.5)*t)
		rz = (1/(2*nbar + 1)) - 2*(p - nbar/(2*nbar + 1)) * np.exp(-2*gamma*(2*nbar + 1)*t)

		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		rho.append(r_list)
	

	return rho


def Hamiltoniano_Sistema(w0):
    
    return (w0/2)*sigmaz()


def Free_Energy(rhot, H, T, St):
    
    rhot_matrix = (1/2)*Qobj([[ 1+rhot[2] , complex(rhot[0],-rhot[1]) ],[ complex(rhot[0],rhot[1]) , 1-rhot[2] ]])
    
    Fneqt = (H*rhot_matrix).tr() - T*St
    
    return Fneqt.real

### MAIN ###

Sr = 0.1

modo = 'Cooling'

Tbanho = 1.442695040888963
w = 2
Tqubit = 4.700782656251331 #0.6606763453727934 #
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))

gamma = 3

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)

curvas, cmodlist = np.loadtxt(f'./FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True)

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

H = Hamiltoniano_Sistema(w0)

Fneq_fatiat = [[],[],[],[],[]]

for i in range(len(cmodlist)):
    
    Fneq = []

    print(cmodlist[i])
    
    clist = np.loadtxt(f'./FisherInformation_{modo}_{Sr}/c_curve_{int(curvas[i])}.txt', unpack=True, dtype=complex, ndmin=1)
    
    rho = RHO(tlist, clist[0], p, gamma, w, nbar)
    
    tlist, S = np.loadtxt(f'./Entropy_{modo}_{Sr}/entropy-{cmodlist[i]:.3f}.txt', unpack=True)    
    
    for t in range(len(rho)):
        
        Fneq.append(Free_Energy(rho[t], H, Tbanho, S[t]))
    
    
    Fneq_fatiat[0].append(Fneq[0])
    Fneq_fatiat[1].append(Fneq[2])
    Fneq_fatiat[2].append(Fneq[5])
    Fneq_fatiat[3].append(Fneq[10])
    Fneq_fatiat[4].append(Fneq[-1])    
        
    
    c = next(colors)
    
    plt.plot(tlist, Fneq, color=c, label=f'|c| = {cmodlist[i]:.3f}')

    
plt.ylabel(r'$F_{neq} (\rho (t))$')
plt.xlabel('Time')
plt.title(f'{modo}')
plt.legend()
plt.xlim(left=0)
plt.tight_layout()
plt.show()


plt.plot(cmodlist, Fneq_fatiat[0], color='black', marker='*', label=f'Time = {tlist[0]:.2f}')
plt.plot(cmodlist, Fneq_fatiat[1], color='black', marker='^', label=f'Time = {tlist[2]:.2f}')
plt.plot(cmodlist, Fneq_fatiat[2], color='black', linestyle=':', label=f'Time = {tlist[5]:.2f}')
plt.plot(cmodlist, Fneq_fatiat[3], color='black', linestyle='--', label=f'Time = {tlist[10]:.2f}')
plt.plot(cmodlist, Fneq_fatiat[4], color='black', linestyle='-', label=f'Time = {tlist[-1]:.2f}')
plt.xlabel('|c|')
plt.ylabel(r'$F_{neq} (\rho)$')
plt.title(f'{modo}')
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()


