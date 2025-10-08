import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import math
import os



def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))

def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(beta*w) + 1)



def RHO(tlist, c, p, gamma, w, nbar, pfinal):
	
	rho = []
	rho_derivada = []
	Srt = []
    
	for t in tlist:

		rx = (c.real/2) * np.exp(-gamma*t)
		ry = -(c.imag/2) * np.exp(-gamma*t)
		rz = (1 - 2*nbar) + 2*(p + nbar - 1) * np.exp(-2*gamma*t)
		
		rmod2 = rx**2 + ry**2 + rz**2

		r_list = [rx, ry, rz, rmod2]
        
		drx = -gamma*(c.real/2) * np.exp(-gamma*t)
		dry = gamma*(c.imag/2) * np.exp(-gamma*t)
		drz = -4 * gamma * (p + nbar - 1) * np.exp(-2*gamma*t)
        
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



def Interseccao_Inicial(metade, pt, pList, Sr_init):
    
    # diferença entre as funções
    f_diff = lambda p: Entropia_Relativa_Populacoes(p, pt) - Sr_init

    # discretiza o intervalo para localizar onde o sinal muda
    Sr_diff = f_diff(np.array(pList))
    
    
    if metade == 0:
        iinit = 0
        ifinal = int(len(pList)/2)
    else:
        iinit = int(len(pList)/2)
        ifinal = len(pList) - 1
        
    print(iinit, ifinal)
    
    for i in range(iinit, ifinal, 1):
          
        if Sr_diff[i] * Sr_diff[i+1] < 0:  # houve cruzamento

            pi = brentq(f_diff, pList[i], pList[i+1])  # raiz exata
            Sri = Entropia_Relativa_Populacoes(pi, pt)
            
    return pi, Sri
    
    
def Interseccao_Temperatura_Inicial(w0, p, Tlist):
    
    find_T = lambda T: pFunc(T, w0) - p
        
    Ti = brentq(find_T, Tlist[0], Tlist[-1])  # busca em todo o range de T
    
    return Ti



def Temperaturas_e_Populacoes(w0, Tlist, p_final, Sr_inicial):
    
    pList = [pFunc(T, w0) for T in Tlist]
    
    plt.plot(Tlist, pList)
    plt.show()
    Sr = []
    
    for i, p in enumerate(pList):
        
        Sr_p = Entropia_Relativa_Populacoes(p, p_final)
        Sr.append(Sr_p)
    
    T_w = Interseccao_Temperatura_Inicial(w0, p_final, Tlist)
    
    p_c, Sr_c = Interseccao_Inicial(0, p_final, pList, Sr_inicial)
    
    T_c = Interseccao_Temperatura_Inicial(w0, p_c, Tlist)
    
    p_h, Sr_h = Interseccao_Inicial(1, p_final, pList, Sr_inicial)
    
    T_h = Interseccao_Temperatura_Inicial(w0, p_h, Tlist)
    
    plt.scatter([p_c], [Sr_c], color='blue')
    plt.scatter([p_h], [Sr_h], color='red')
    plt.plot(pList, Sr, label=f'pf = {p_final}', color='orange')
    plt.hlines(y=Sr_inicial, xmin=min(pList), xmax=max(pList), color='black')
    plt.xlabel('Populations')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()

    return p_c, p_h, Sr_c, Sr_h, Tc, Th, Tw



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

p_final = 0.8

Sr_inicial = 0.1 ### se alterar, deve mudar os ranges de temperatura em Temperaturas_e_Populacoes

tlist = np.arange(0, 10, 0.01)

print('Temperatura e Populações')

TList = np.arange(0.1, 100, 0.0001)

p_c, p_h, Sr_c, Sr_h, Tc, Th, Tw = Temperaturas_e_Populacoes(w0, TList, p_final, Sr_inicial)

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














