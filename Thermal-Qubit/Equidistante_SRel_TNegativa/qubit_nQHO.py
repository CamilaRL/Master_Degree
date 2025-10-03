import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import math
import os



def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) - 1)



def RHO(tlist, c, p, gamma, w, nbar, pfinal):
	
	rho = []
	rho_derivada = []
	Srt = []
    
	for t in tlist:

		rx = c.real * np.exp(-2*gamma*(nbar + 0.5)*t)
		ry = -c.imag * np.exp(-2*gamma*(nbar + 0.5)*t)
		rz = -(nbar/(2*nbar + 1)) - (p - (nbar + 1)/(2*nbar + 1)) * np.exp(-2*gamma*(2*nbar + 1)*t)
		
		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		drx = c.real * (-2*gamma*(nbar + 0.5)) * np.exp(-2*gamma*(nbar + 0.5)*t)
		dry = -c.imag * (-2*gamma*(nbar + 0.5)) * np.exp(-2*gamma*(nbar + 0.5)*t)
		drz = 2 * gamma * (2*nbar + 1) *(p - (nbar + 1)/(2*nbar + 1)) * np.exp(-2*gamma*(2*nbar + 1)*t)
        
		drmod2 = drx**2 + dry**2 + drz**2
        
		dr_list = [drx, dry, drz, drmod2]
        
		rho.append(r_list)
		rho_derivada.append(dr_list)
        
		Srt.append(Entropia_Relativa_Bloch(rho[0], rho[-1]))

	return rho, rho_derivada, Srt

def Entropia_Relativa_BP(rho_i, p_final):

    autoval_i = [(1 + np.sqrt(rho_i[3]))/2, (1 - np.sqrt(rho_i[3]))/2]
    autoval_f = [p_final, 1 - p_final]
    
    Sr = 0
    
    for k in range(2):
        Sr = Sr + autoval_i[k] * np.log(autoval_i[k] / autoval_f[k])
    
    return Sr


def Entropia_Relativa_Bloch(rho_i, rho_f):

    autoval_i = [(1 + np.sqrt(rho_i[3]))/2, (1 - np.sqrt(rho_i[3]))/2]
    autoval_f = [(1 + np.sqrt(rho_f[3]))/2, (1 - np.sqrt(rho_f[3]))/2]
    print(np.sqrt(rho_f[3]))
    Sr = 0
    
    for k in range(2):
        Sr = Sr + autoval_i[k] * np.log(autoval_i[k] / autoval_f[k])
    
    return Sr


def Entropia_Relativa_Populacoes(pi, pt):

    return pi*np.log(pi/pt) + (1-pi)*np.log((1-pi)/(1-pt))
    


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



def Intersection_Initial_Sr(pt, pList, Tlist, w0, Sr_init):
    
    # diferença entre as funções
    f_diff = lambda p: Entropia_Relativa_Populacoes(p, pt) - Sr_init

    # discretiza o intervalo para localizar onde o sinal muda
    y_diff = f_diff(np.array(pList))

    xs, ys = [], []
    Ts = []
    
    for i in range(len(pList) - 1):
        
        if y_diff[i] * y_diff[i+1] < 0:  # houve cruzamento
            
            xi = brentq(f_diff, pList[i], pList[i+1])  # raiz exata
            yi = Entropia_Relativa_Populacoes(xi, pt)
            
            xs.append(xi)
            ys.append(yi)
            
            # agora inverte para achar T correspondente
            g = lambda T: pFunc(T, w0) - xi
            T_root = brentq(g, Tlist[0], Tlist[-1])  # busca em todo o range de T
            Ts.append(T_root)
            
    g = lambda T: pFunc(T, w0) - pt
    Tw = brentq(g, Tlist[0], Tlist[-1])

    return xs, ys, Ts[0], Ts[1], Tw



def Temperaturas_e_Populacoes(w0, p_final, Sr_inicial):
    
    Tlist = np.arange(-10, 1, 0.00001)
    
    pList = [pFunc(T, w0) for T in Tlist]
    
    Sr = []

    for i, p in enumerate(pList):
        
        Sr_p = Entropia_Relativa_Populacoes(p, p_final)
        Sr.append(Sr_p)
    
    
    pinit, Sinit, Tc, Th, Tw = Intersection_Initial_Sr(p_final, pList, Tlist, w0, Sr_inicial)
    
    plt.scatter(pinit, Sinit)
    plt.plot(pList, Sr, label=f'pf = {p_final}')
    plt.hlines(y=Sr_inicial, xmin=min(pList), xmax=max(pList), color='black')
    plt.xlabel('Populations')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()

    return pinit[0], pinit[1], Sinit, Tc, Th, Tw



def WriteOutput(processo, tlist, FQ_classificados, Sr):
    
    os.mkdir(f'./FisherInformation_{processo}')

    f_cmod = open(f'./FisherInformation_{processo}/cmod.txt', 'w')

    f_cmod.write(f'0 0\n')
    
    f_cmod.close()

    i = 0
    
    f_fisher = open(f'./FisherInformation_{processo}/QFI_curve_{i}.txt', 'w')
    
    f_Sr = open(f'./FisherInformation_{processo}/Sr_curve_{i}.txt', 'w')

    for t in range(len(tlist)):
        
        f_fisher.write(f'{tlist[t]} {FQ_classificados[t]}\n')
            
        f_Sr.write(f'{tlist[t]} {Sr[t]}\n')
    
    f_fisher.close()
    f_Sr.close()


##### MAIN #####

## parametros

w = 2
w0 = 2

gamma = 3

p_final = 0.5

Sr_inicial = 0.1

tlist = np.arange(0, 10, 0.01)

print('Temperatura e Populações')
pc, ph, Sinit, Tc, Th, Tw = Temperaturas_e_Populacoes(w0, p_final, Sr_inicial)

nbar = nbarFunc(Tw, w)

### prints de checagem

print(f'pc {Tc} {pc} - {pFunc(Tc, w0)}')
print(f'ph {Th} {ph} - {pFunc(Th, w0)}')
print(f'pw {Tw} {p_final} - {pFunc(Tw, w0)}')

print(f'Sc {Sinit[0]} - {Entropia_Relativa_Populacoes(pFunc(Tc, w0), p_final)}')
print(f'Sh {Sinit[1]} - {Entropia_Relativa_Populacoes(pFunc(Th, w0), p_final)}')

## Aquecer

print('Aquecer')

rhot, drhot, Srt_aquecer = RHO(tlist, 0, pc, gamma, w, nbar, p_final)
    
QFI_aquecer = FisherInformation(rhot, drhot)


## Resfriamento

print('Resfriar')
    
rhot, drhot, Srt_resfriar = RHO(tlist, 0, ph, gamma, w, nbar, p_final)
    
QFI_resfriar = FisherInformation(rhot, drhot)
    



## write files ##

WriteOutput(f'Heating_{Sr_inicial}', tlist, QFI_aquecer, Srt_aquecer)
WriteOutput(f'Cooling_{Sr_inicial}', tlist, QFI_resfriar, Srt_resfriar)














