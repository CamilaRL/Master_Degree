import numpy as np
import matplotlib.pyplot as plt
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



def Temperatura_Efetiva(rho, w0):
    
    lamb_mais = 1 + np.sqrt(rho[3])
    lamb_menos = 1 - np.sqrt(rho[3])
    
    E_mais = -w0/2
    E_menos = w0/2
    
    Teff = (E_menos - E_mais)/np.log(lamb_mais/lamb_menos)

    return Teff
    
    
### MAIN ###

Tbanho = 2
w = 2
Tqubit = 10
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
    
    Teff = []
    for t in range(len(rho)):
        Teff.append(Temperatura_Efetiva(rho[t], w0))
    
    c = next(colors)
    plt.plot(tlist, Teff, color=c, label=f'|c| = {cmodlist[i]:.3f}')
    
plt.hlines(y=Tbanho, xmin=tlist[0], xmax=tlist[-1], linestyle='--', color='black', label='Thermal Bath')
plt.hlines(y=Tqubit, xmin=tlist[0], xmax=tlist[-1], linestyle=':', color='black', label='Qubit')
plt.ylabel('Temperature')
plt.xlabel('Time')
plt.legend(loc='center right')
plt.tight_layout()
plt.show()