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
		rz = -(1/(2*nbar + 1)) + 2*(p - nbar/(2*nbar + 1)) * np.exp(-2*gamma*(2*nbar + 1)*t)

		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		rho.append(r_list)
	

	return rho


def Temperatura_Real(rho, w0):
    
    pg = (1 + rho[2])/2
    pe = (1 - rho[2])/2

    T = w0/np.log(pg/pe)

    return T

### MAIN ###

modo = 'Aquecer'

Tbanho = 10
w = 2
Tqubit = 2
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))

gamma = 3

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)

curvas, cmodlist = np.loadtxt(f'./FisherInformation_{modo}/cmod.txt', unpack=True)

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))


for i in range(len(cmodlist)):
    
    Treal = []

    print(cmodlist[i])
    
    clist = np.loadtxt(f'./FisherInformation_{modo}/c_curve_{int(curvas[i])}.txt', unpack=True, dtype=complex, ndmin=1)
    
    rho = RHO(tlist, clist[0], p, gamma, w, nbar)
    
    for t in range(len(tlist)):
        Treal.append(Temperatura_Real(rho[t], w0))
    
    
    c = next(colors)
    
    plt.plot(tlist, Treal, color=c, label=f'|c| = {cmodlist[i]:.3f}')


plt.hlines(y=Tbanho, xmin=tlist[0], xmax=tlist[-1], linestyle='--', color='black', label='Thermal Bath')
plt.hlines(y=Tqubit, xmin=tlist[0], xmax=tlist[-1], linestyle=':', color='black', label='Qubit')
plt.ylabel('Real Temperature')
plt.xlabel('Time')
plt.title('Heating')
plt.legend(loc='upper right', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()

    
    
    
    
    
    
