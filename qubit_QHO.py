import numpy as np
import matplotlib.pyplot as plt
from qutip import *



def DERIVADA_RHO(tlist, c, p, g, w, nbar):

	rho = []
	
	for t in tlist:
	
		coerencia = c*4*g*g.conjugate()*(2*nbar + 1)*np.sin(w*t)*np.exp(4*g*g.conjugate()*(2*nbar + 1)*(np.cos(w*t) - 1)/w**2)/w
	
		rhot = Qobj( [[0, -coerencia], [-coerencia.conjugate(), 0]] )
		
		rho.append(rhot)
	
	return rho


def FisherInformation(rho, drho):

	
	FQ = []
	FQimg = []

	for t in range(len(rho)):

		aval, avet = rho[t].eigenstates()
		
		FQt = 0
		
		for i in range(len(aval)):
	
			for j in range(len(aval)):
		
				den = aval[i] + aval[j]
				
				if abs(den) > 1e-12:
				
					sanduiche = (avet[i].dag() * drho[t] * avet[j])[0,0]

					FQt = FQt + (2*sanduiche*(sanduiche.conjugate()) / (den))
		
		FQ.append(FQt.real)
		FQimg.append(FQt.imag)
		
	return FQ, FQimg
		

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(-beta*w) - 1)

def RHO(tlist, c, p, g, w, nbar):

	rho = []
	
	for t in tlist:
        
		coerencia = Verifica_Positividade(t, c, p, g, w, nbar)
	
		rhot = Qobj( [[p, coerencia], [coerencia.conjugate(), 1-p]] )
        
		rho.append(rhot)
	
	return rho


def Verifica_Positividade(t, c, p, g, w, nbar):
    
    coerencia = c*np.exp(4*g*g.conjugate()*(2*nbar + 1)*(np.cos(w*t) - 1)/w**2)
    
    coerencia_max = np.sqrt(max(0.0, p*(1-p)))
    
    mag = abs(coerencia)
    
    if mag > coerencia_max:
        
        coerencia = (coerencia/(mag + 0j)) * coerencia_max
    
    return coerencia


### MAIN ###

T = 10
w = 2

c = complex(0.5, 0.9)
p = 0.5
g = complex(0.2, 0.1)

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(T, w)

rho = RHO(tlist, c, p, g, w, nbar)

drho = DERIVADA_RHO(tlist, c, p, g, w, nbar)

FQlist, FQlistIMG = FisherInformation(rho, drho)


plt.plot(tlist, FQlist)
plt.ylabel('Fisher Information')
plt.xlabel('Time')
plt.tight_layout()
plt.show()






















