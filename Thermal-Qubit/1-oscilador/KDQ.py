import numpy as np
import matplotlib.pyplot as plt

def RHO(tlist, c, p, g, w, nbar):
	
	rho = []
	rho_derivada = []
	
	for t in tlist:
		
		coerencia = c*np.exp(4*g*g.conjugate()*(2*nbar + 1)*(np.cos(w*t) - 1)/w**2)
		
		fator = -4*g*g.conjugate()*(2*nbar + 1)*np.sin(w*t)/w
		coerencia_derivada = coerencia * fator
		
		rhot = Qobj( [[p, coerencia], [coerencia.conjugate(), 1-p]] )
		
		rhot_derivada = Qobj( [[0, coerencia_derivada], [coerencia_derivada.conjugate(), 0]] )
        
		rho.append(rhot)
		rho_derivada.append(rhot_derivada)
	

	return rho, rho_derivada
	
	
def Coerencia_Positividade(tlist, c, p, g, w, nbar):
    
    coerencia_max = np.sqrt(max(0.0, p*(1-p)))   
     
    clist = []
    
    for t in tlist:
    	
    	cnew = c
    	
    	coerencia = cnew*np.exp(4*g*g.conjugate()*(2*nbar + 1)*(np.cos(w*t) - 1)/w**2)
       
    	mag = abs(coerencia)
    
    	if mag > coerencia_max:
        
        	cnew = (cnew/(mag + 0j)) * coerencia_max
       	
    	clist.append(cnew)
    
    clist_mod = []
    for ci in clist:
        clist_mod.append(abs(ci))
    
    cfinal = clist[clist_mod.index(min(clist_mod))]

    return cfinal
