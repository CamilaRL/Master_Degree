import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
import math
import os



def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))

    
def Entropia_Relativa_Populacoes(pi, pt):
    
    return pi*np.log(pi/pt) + (1-pi)*np.log((1-pi)/(1-pt))


def Interseccao_Inicial(metade, pt, pList, Sr_init):
    
    # diferença entre as funções
    f_diff = lambda p: Entropia_Relativa_Populacoes(p, pt) - Sr_init

    # discretiza o intervalo para localizar onde o sinal muda
    Sr_diff = f_diff(np.array(pList))
    
    if metade == 0:
        iinit = 0
        ifinal = np.argmin(Sr_diff)
    else:
        iinit = np.argmin(Sr_diff)
        ifinal = len(pList) - 1
        
    for i in range(iinit, ifinal, 1):

        if Sr_diff[i] * Sr_diff[i+1] < 0:  # houve cruzamento

            pi = brentq(f_diff, pList[i], pList[i+1])  # raiz exata
            Sri = Entropia_Relativa_Populacoes(pi, pt)
            
    return pi, Sri
    
    
def Interseccao_Temperatura_Inicial(w0, p, Tlist):
    
    find_T = lambda T: pFunc(T, w0) - p
        
    Ti = brentq(find_T, Tlist[0], Tlist[-1])  # busca em todo o range de T
    
    return Ti


def Temperaturas_e_Populacoes(metade, w0, Tlist, p_final, Sr_inicial):
    
    pList = [pFunc(T, w0) for T in Tlist]
    
    Sr = []
    
    for i, p in enumerate(pList):
        
        Sr_p = Entropia_Relativa_Populacoes(p, p_final)
        Sr.append(Sr_p)
    
    p_i, Sr_i = Interseccao_Inicial(metade, p_final, pList, Sr_inicial)
    
    T_i = Interseccao_Temperatura_Inicial(w0, p_i, Tlist)

    return p_i, Sr_i, T_i, Sr, pList


def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f


def Coherences(w0, beta_1, beta_2, Z1, Z2):

    alpha_list = []
    mod_list = []

    alpha = np.exp(-w0*(beta_1 + beta_2)/2)/(Z1*Z2)
    
    num = np.arange(0.01, 1, 0.01)

    for r in num:
    
        for i in num:
        
            c = complex(r, i)
            c_abs = abs(c)
            
            if (c_abs <= alpha) and (c_abs not in mod_list):
                
                alpha_list.append(c)

                mod_list.append(c_abs)
    
    return alpha_list


def Master_Equation(w0, beta_1, beta_2, tlist, L_operators):

    ## density matrices

    Z1 = 2*np.cosh(beta_1*w0)
    Z2 = 2*np.cosh(beta_2*w0)

    p1 = np.exp(beta_1*w0)/Z1
    p2 = np.exp(beta_2*w0)/Z2

    rho0_q1 = Qobj([[p1, 0],[0, 1-p1]])

    rho0_q2 = Qobj([[p2, 0],[0, 1-p2]])


    alpha_list = Coherences(w0, beta_1, beta_2, Z1, Z2)

    G = basis(2,0)
    E = basis(2,1)

    v00 = tensor(G, G)
    v01 = tensor(G, E)
    v10 = tensor(E, G)
    v11 = tensor(E, E)

    #coherences_matrix = alpha_list[0] * v01 * v10.dag() + alpha_list[0].conjugate() * v10 * v01.dag()

    rho0 = tensor(rho0_q1, rho0_q2) #+ coherences_matrix


    ## solve master equation

    dataME = mesolve(H_S, rho0, tlist, L_operators, [])

    rhof = dataME.states
    
    rhof_q1 = []
    rhof_q2 = []
    
    for rt in rhof:
    
        rhof_q1.append( rt.ptrace([0]) )
        rhof_q2.append( rt.ptrace([1]) )

    return rhof, rhof_q1, rhof_q2


def Write_Density_Matrices(rhof, rhof_q1, rhof_q2, modo):

    os.mkdir(f'./DensityMatrices_c0_{modo}')

    ## write rho in output filess

    for t, rhot in enumerate(rhof):
        
        rhot_q1 = rhof_q1[t]
        rhot_q2 = rhof_q2[t]

        f = open(f'./DensityMatrices_c0_{modo}/rhof_t{t}.txt', 'w')
        f1 = open(f'./DensityMatrices_c0_{modo}/rhof_q1_t{t}.txt', 'w')
        f2 = open(f'./DensityMatrices_c0_{modo}/rhof_q2_t{t}.txt', 'w')

        for i in range(4):
            for j in range(4):
            
                f.write(f'{rhot[i][j]} ')
                
            f.write('\n')
        
        for i in range(2):
            for j in range(2):
            
                f1.write(f'{rhot_q1[i][j]} ')
                f2.write(f'{rhot_q2[i][j]} ')

            f1.write('\n')
            f2.write('\n')
                
    f.close()
    f1.close()
    f2.close()

    


######################## MAIN ########################

## parameters

w0 = 2
gamma = 1

tlist = np.arange(0.005, 2, 0.001)

p_final = 0.8

Sr_inicial = 0.1


## hamiltonians

H_q1 = w0*sigmaz()/2

H_q2 = w0*sigmaz()/2

H_S = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)


## temperaturas equidistantes

Tc_list = np.arange(0.1, 100, 0.0001)
Th_list = np.arange(0.1, 100, 0.0001)

Tw = Interseccao_Temperatura_Inicial(w0, p_final, Th_list)

pc, Src, Tc, Sr0, pList0 = Temperaturas_e_Populacoes(0, w0, Tc_list, p_final, Sr_inicial)

ph, Srh, Th, Sr1, pList1 = Temperaturas_e_Populacoes(1, w0, Th_list, p_final, Sr_inicial)

print(Tc, Tw, Th)

plt.plot(pList0, Sr0, color='orange', label=f'pf = {p_final}')
plt.plot(pList1, Sr1, color='orange')
plt.scatter([p_final], [0], color='orange', label=f'Tw = {Tw}')
plt.scatter([pc], [Src], color='blue', label=f'Tc = {Tc}')
plt.scatter([ph], [Srh], color='red', label=f'Th = {Th}')
plt.hlines(y=Sr_inicial, xmin=min(pList1), xmax=max(pList0), color='black', label='Initial Relative Entropy')
plt.xlabel('Populations')
plt.ylabel('Relative Entropy')
plt.legend()
plt.show()


beta_c = 1/Tc
beta_h = 1/Th
beta_R = 1/Tw


## collapse operators

L_operators = []

evals, evecs = H_S.eigenstates()

for n, valn in enumerate(reversed(evals)):
    for m, valm in enumerate(reversed(evals)):
        
        if m < n:
            
            dmn = valm - valn
            
            if dmn != 0:
                
                f = Distribution('bose', beta_R, dmn)
                
                C_abs = gamma * np.sqrt(1 + f) * (evecs[n] * evecs[m].dag())
                
                L_operators.append(C_abs)
                
                C_emi = gamma * np.sqrt(f) * (evecs[m] * evecs[n].dag())

                L_operators.append(C_emi)


## master equation solver

rhof_heating, rhof_heating_q1, rhof_heating_q2 = Master_Equation(w0, beta_c, beta_c, tlist, L_operators)

rhof_cooling, rhof_cooling_q1, rhof_cooling_q2 = Master_Equation(w0, beta_h, beta_h, tlist, L_operators)


## write output file

Write_Density_Matrices(rhof_heating, rhof_heating_q1, rhof_heating_q2, 'Heating')

Write_Density_Matrices(rhof_cooling, rhof_cooling_q1, rhof_cooling_q2, 'Cooling')










