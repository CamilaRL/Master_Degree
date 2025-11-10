import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math as m



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
				
					sanduiche = (avet[i].dag() * drho[t] * avet[j])

					FQt = FQt + (2*sanduiche*(sanduiche.conjugate()) / den)
		
		FQ.append(FQt.real)
		
		if FQt.imag != 0:
		    print('Informação de Fisher Imaginária!')

		
	return FQ
		

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) - 1)


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


def Classifica_FQ(FQ_classificados, FQ, classes, c):

    inicio_curva = [FQ[1], FQ[2], FQ[3], FQ[4], FQ[5]]
    
    inicio_curva_class = []
    
    for curva in FQ_classificados:
        
        inicio_curva_class.append([curva[1], curva[2], curva[3], curva[4], curva[5]])
    
    mesma_curva = False
    
    
    for idx, curva in enumerate(FQ_classificados):
    
        inicio_curva_class = [curva[1], curva[2], curva[3], curva[4], curva[5]]
    
        if all(m.isclose(inicio_curva[i], inicio_curva_class[i]) for i in range(5)):

            mesma_curva = True
            k = idx
            break
                
    if mesma_curva:
        
        classes[k].append(c)

    else:
        
        FQ_classificados.append(FQ)
        classes.append([c])
        

    return FQ_classificados, classes


def FileWrite(path, QFI, tlist, classe, curva):

    f = open(path+f'curva_{curva}.txt', 'w')
    cfile = open(path+f'c_curva_{curva}.txt', 'w')
    
    for i in range(len(tlist)):
        
        f.write(f'{tlist[i]} {QFI[i]}\n')
    
    for c in classe:
    
        cfile.write(f'{c}\n')
        
    f.close()
    cfile.close()
    


### MAIN ###

Tbanho = 2
w = 2
Tqubit = 10
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))
g = complex(0.1, 0.1)

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)


clist = []
clist_mod = []
ilist = np.arange(0, 1, 0.1)
for i in ilist:
    for j in ilist:
    
        c = complex(i, j)
        cnew = Coerencia_Positividade(tlist, c, p, g, w, nbar)
        
        if cnew not in clist:
            clist.append(cnew)
            
clist = np.sort_complex(clist)


FQ_classificados = []
classes = []
maximos = []

for c in clist:
    
    print(c)
    
    rho, drho = RHO(tlist, c, p, g, w, nbar)
    
    FQ = FisherInformation(rho, drho)
    
    FQ_classificados, classes = Classifica_FQ(FQ_classificados, FQ, classes, c)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, len(classes))))

for i in range(len(classes)):

    cor = next(colors)
    
    nome = f'Curve {i}'
    
    print(f'{nome}')
    
    FileWrite('./FisherInformation/', FQ_classificados[i], tlist, classes[i], i)
    
    plt.scatter(tlist, FQ_classificados[i], color=cor, s=1, label=f'{nome}')


plt.ylabel('Fisher Information')
plt.xlabel('Time')
plt.legend()
plt.tight_layout()
plt.show()




















