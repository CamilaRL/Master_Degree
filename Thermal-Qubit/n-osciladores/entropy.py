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

Tbanho = 10
w = 2
Tqubit = 2
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

Slist = []

for i in range(len(cmodlist)):
    
    print(cmodlist[i])
    
    rho = RHO(tlist, clist[cindex[i][0]], p, gamma, w, nbar)
    
    S = vonNeumann_Entropy(rho)
    
    equilibrio = np.where(S==S[-1])
    
    print(f'{clist[cindex[i][0]]:.3f} : t={tlist[equilibrio[0][0]]} ii={equilibrio[0][0]} if={equilibrio[0][-1]}')
    
    Slist.append(S)
    
    c = next(colors)
    
    #WriteFile(f'./Entropy/entropy-{cmodlist[i]:.3f}.txt', S, tlist)
    
    plt.plot(tlist, S, color=c, label=f'|c| = {cmodlist[i]:.3f}')

plt.ylabel('S')
plt.xlabel('Time')
plt.title('Cooling')
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
plt.title('Heating')
plt.xscale('log')
plt.xlim(left=0.01)
plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()
