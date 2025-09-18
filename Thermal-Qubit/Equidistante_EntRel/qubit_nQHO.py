import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math




def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) - 1)


def Coerencia(p):

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
    
    return clist


def RHO(tlist, c, p, gamma, w, nbar):
	
	rho = []
	rho_derivada = []
	
	for t in tlist:

		rx = c.real * np.exp(-2*gamma*(nbar + 0.5)*t)
		ry = -c.imag * np.exp(-2*gamma*(nbar + 0.5)*t)
		rz = (1/(2*nbar + 1)) - 2*((nbar + 1)/(2*nbar + 1) - p) * np.exp(-2*gamma*(2*nbar + 1)*t)

		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		drx = c.real * (-2*gamma*(nbar + 0.5)) * np.exp(-2*gamma*(nbar + 0.5)*t)
		dry = -c.imag * (-2*gamma*(nbar + 0.5)) * np.exp(-2*gamma*(nbar + 0.5)*t)
		drz = 4 * gamma * (2*nbar + 1) * ((nbar + 1)/(2*nbar + 1) - p) * np.exp(-2*gamma*(2*nbar + 1)*t)
        
		drmod2 = drx**2 + dry**2 + drz**2
        
		dr_list = [drx, dry, drz, drmod2]
        
		rho.append(r_list)
		rho_derivada.append(dr_list)
	

	return rho, rho_derivada


def Entropia_Relativa(rho_i, rho_f):

    autoval_i = [1 + np.sqrt(rho_i[3]), 1 - np.sqrt(rho_i[3])]
    autoval_f = [1 + np.sqrt(rho_f[3]), 1 - np.sqrt(rho_f[3])]
    
    Sr = 0
    
    for k in range(2):
        Sr = Sr + autoval_i[k] * np.log(autoval_i[k] / autoval_f[k])
    
    return Sr


def FisherInformation(rho, drho):
	
	FQ = []

	for t in range(len(rho)):

		FQt = drho[t][3] + ((rho[t][0]*drho[t][0] + rho[t][1]*drho[t][1] + rho[t][2]*drho[t][2])**2)/(1 - rho[t][3])
		
		FQ.append(FQt)
        
		if FQt.imag != 0:
		    print('Informação de Fisher Imaginária!')
		elif abs(FQt) < 0:
			print('Informação de Fisher Negativa!')
		
	return FQ


def Classifica_FQ(FQ_classificados, FQ, c_classes_list, cmod_list, c):

    cmod = abs(c)
    
    novo = True
    
    for i, cmodi in enumerate(cmod_list):
    
        if math.isclose(cmodi, cmod, abs_tol=0.00001):
            novo = False
            c_classes_list[i].append(c)
    
    if novo:
        
        FQ_classificados.append(FQ)
        
        cmod_list.append(cmod)
        
        c_classes_list.append([c])
        

    return FQ_classificados, c_classes_list, cmod_list



##### MAIN #####

## parametros

Th = 500
Tw = 5
Tc = 2

w = 2
w0 = 2

gamma = 3

nbar = nbarFunc(Tw, w)

tlist = np.arange(0, 10, 0.01)

pc = np.exp(w0/(2*Tc))/(2*np.cosh(w0/(2*Tc)))

clist_c = Coerencia(pc)


## Aquecer

rho_aquecer = []
drho_aquecer = []
Sr_aquecer = []

for c in clist_c:
    
    rho, drho = RHO(tlist, c, pc, gamma, w, nbar)
    
    rho_aquecer.append(rho)
    drho_aquecer.append(drho)
    
    Sr_aquecer.append(Entropia_Relativa(rho[0], rho[-1]))
    

## variacao de Th para ter mesma entropia relativa que Tc 
Sr = 0
while(Sr_aquecer[0] != Sr):

    Th = Th + 10
    
    ph = np.exp(w0/(2*Th))/(2*np.cosh(w0/(2*Th)))
    
    rho, drho = RHO(tlist, 0, ph, gamma, w, nbar)

    Sr = Entropia_Relativa(rho[0], rho[-1])
    
    print(f'{Sr_aquecer[0]} {Th} {Sr}')

'''
## Resfriamento

rho_resfriar = []
drho_resfriar = []
Sr_resfriar = []

for c in clist_h:
    
    rho, drho = RHO(tlist, c, ph, gamma, w, nbar)
    
    rho_resfriar.append(rho)
    drho_resfriar.append(drho)
    
    Sr_resfriar.append(Entropia_Relativa(rho[0], rho[-1]))
    


'''
        



'''
## write files ##


curvas = np.arange(0, len(cmod_list), 1, dtype=np.int16)

cmod_list, curvas = (list(t) for t in zip(*sorted(zip(cmod_list, curvas))))


f_cmod = open(f'./FisherInformation/cmod.txt', 'w')

for i in range(len(curvas)):
    
    f_cmod.write(f'{curvas[i]} {cmod_list[i]}\n')
    
f_cmod.close()


for i in curvas:
    
    f_fisher = open(f'./FisherInformation/QFI_curve_{i}.txt', 'w')
    
    f_c_classes = open(f'./FisherInformation/c_curve_{i}.txt', 'w')


    for t in range(len(tlist)):
    
        f_fisher.write(f'{tlist[t]} {FQ_classificados[i][t]}\n')
    
    for c in c_classes[i]:
    
        f_c_classes.write(f'{c}\n')
        
    f_fisher.close()
    f_c_classes.close()


'''
















