import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math



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


def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) - 1)


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


FQ_classificados = []
c_classes = []
cmod_list = []

for c in clist:
    
    print(c)
    
    rho, drho = RHO(tlist, c, p, gamma, w, nbar)
    
    FQ = FisherInformation(rho, drho)
    
    FQ_classificados, c_classes, cmod_list = Classifica_FQ(FQ_classificados, FQ, c_classes, cmod_list, c)


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



















