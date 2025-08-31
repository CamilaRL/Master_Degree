import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math as m



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
		rz = (1/(2*nbar + 1)) - 2*(p - nbar/(2*nbar + 1)) * np.exp(-2*gamma*(2*nbar + 1)*t)

		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		drx = c.real * (-2*gamma*(nbar + 0.5)) * np.exp(-2*gamma*(nbar + 0.5)*t)
		dry = -c.imag * (-2*gamma*(nbar + 0.5)) * np.exp(-2*gamma*(nbar + 0.5)*t)
		drz = 2 * (2*gamma*(2*nbar + 1)) * (p - nbar/(2*nbar + 1)) * np.exp(-2*gamma*(2*nbar + 1)*t)
        
		drmod2 = drx**2 + dry**2 + drz**2
        
		dr_list = [drx, dry, drz, drmod2]
        
		rho.append(r_list)
		rho_derivada.append(dr_list)
	

	return rho, rho_derivada


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

Tbanho = 10
w = 2
Tqubit = 19
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
classes = []
maximos = []

for c in clist:
    
    print(c)
    
    rho, drho = RHO(tlist, c, p, gamma, w, nbar)
    
    FQ = FisherInformation(rho, drho)
    
    FQ_classificados, classes = Classifica_FQ(FQ_classificados, FQ, classes, c)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.1, 1, len(classes))))

for i in range(len(classes)):

    cor = next(colors)
    
    nome = f'Curve {i}'
    
    print(f'{nome}')
    
    FileWrite('./FisherInformation/', FQ_classificados[i], tlist, classes[i], i)
    
    plt.scatter(tlist, FQ_classificados[i], color=cor, s=10, label=f'{nome}')


plt.ylabel('Fisher Information')
plt.xlabel('Time')
plt.legend()
plt.tight_layout()
plt.show()


















