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
    
    return Fneqt

### MAIN ###

Tbanho = 10
w = 2
Tqubit = 9
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))

gamma = 3

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)

coerencia_max = np.sqrt(max(0.0, p*(1-p)))

clist = []
ilist = np.arange(0, 1, 0.1)
for i in ilist:
    for j in ilist:
    
        c = complex(i, j)
        
        mag = abs(c)
        
        if mag > coerencia_max:
            c = (c/mag) * coerencia_max
        
        if c not in clist:
            clist.append(c)
            
clist = np.sort_complex(clist)


## Separa c pelo modulo pois curvas com o mesmo modulo tem o mesmo comportamento
cmodlist = [abs(clist[0])]
cindex = [[0]]
for c in range(1, len(clist), 1):
    
    cmod = abs(clist[c])
    
    novo = True
    for i, cmodi in enumerate(cmodlist):
        
        if math.isclose(cmod, cmodi, rel_tol=1e-10):
            novo = False
            cindex[i].append(c)
    
        elif novo and i == (len(cmodlist)-1):
            
            cmodlist.append(cmod)
            cindex.append([c])


cmodlist, cindex = (list(t) for t in zip(*sorted(zip(cmodlist, cindex))))

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

H = Hamiltoniano_Sistema(w0)

for i in range(len(cmodlist)):
    
    Fneq = []

    print(cmodlist[i])
    
    rho = RHO(tlist, clist[cindex[i][0]], p, gamma, w, nbar)
    
    tlist, S = np.loadtxt(f'./Entropy/dT-{abs(Tbanho-Tqubit)}-Aquecer/entropy-{cmodlist[i]:.3f}.txt', unpack=True)    
    
    for t in range(len(rho)):
        
        Fneq.append(Free_Energy(rho[t], H, Tbanho, S[t]))
    
    c = next(colors)
    
    plt.plot(tlist, Fneq, color=c, label=f'|c| = {cmodlist[i]:.3f}')
    
plt.ylabel(r'$F_{neq} (\rho (t))$')
plt.xlabel('Time')
plt.title(r'Heating ($\Delta$T = ' + f'{abs(Tbanho-Tqubit)})')
plt.legend()
plt.xlim(left=0, right=0.15)
plt.tight_layout()
plt.show()