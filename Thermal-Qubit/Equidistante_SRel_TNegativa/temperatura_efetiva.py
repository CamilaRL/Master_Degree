import numpy as np
import matplotlib.pyplot as plt
import math

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) + 1)
    

def RHO(tlist, c, p, gamma, w, nbar):
	
	rho = []
	rho_derivada = []
	
	for t in tlist:

		rx = (c.real/2) * np.exp(-gamma*t)
		ry = -(c.imag/2) * np.exp(-gamma*t)
		rz = (1 - 2*nbar) + 2*(p + nbar - 1) * np.exp(-2*gamma*t)

		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		rho.append(r_list)
	

	return rho



def Temperatura_Efetiva(rho, w0):
    
    lamb_mais = (1 + np.sqrt(rho[3]))/2
    lamb_menos = (1 - np.sqrt(rho[3]))/2
    
    E_mais = -w0/2
    E_menos = w0/2
    
    Teff = (E_menos - E_mais)/np.log(lamb_mais/lamb_menos)

    return Teff
    
    
### MAIN ###

modo = 'Cooling' ################### MUDAR
Sr = 0.1

Tbanho = 4.932606924752464
w = 2
Tqubit = -4.013736266992159 ################### MUDAR
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))

gamma = 1

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)

curvas = np.loadtxt(f'./FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
cmodlist = np.loadtxt(f'./FisherInformation_{modo}_{Sr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)

cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmodlist))))

Teff_inicial = []

for i, cmod in enumerate(cmodlist):
    
    print(cmod)
    
    Teff = []
    
    rho = RHO(tlist, 0, p, gamma, w, nbar)
        
    for t in range(len(rho)):
        Teff.append(Temperatura_Efetiva(rho[t], w0))
    
    Teff_inicial.append(Teff[0])
    
    c = next(colors)
    plt.plot(tlist, Teff, color=c, label=f'|c| = {cmod:.3f}')
    
plt.hlines(y=Tbanho, xmin=tlist[0], xmax=tlist[-1], linestyle='--', color='black', label='Thermal Bath')
plt.hlines(y=Tqubit, xmin=tlist[0], xmax=tlist[-1], linestyle=':', color='black', label='Qubit')
plt.ylabel('Effective Temperature')
plt.xlabel('Time')
plt.title(modo)
plt.legend(loc='upper right', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()

plt.scatter(cmodlist, Teff_inicial, s=15, color='black')
plt.plot(cmodlist, Teff_inicial, color='black')
plt.ylabel('Initial Effective Temperature')
plt.xlabel('|c|')
plt.title(modo)
plt.show()







