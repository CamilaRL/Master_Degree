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
    
    Ti = w0 / np.log(p/(1-p))
    
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


def Coherences(w0, beta_1, beta_2):

    alpha_list = []
    mod_list = []

    Z1 = 2*np.cosh(beta_1*w0/2)
    Z2 = 2*np.cosh(beta_2*w0/2)

    alpha = np.exp(-w0*(beta_1 + beta_2)/2)/(Z1*Z2)
    
    num = np.linspace(0, alpha, 3)
    
    for i in num:
        for j in num:
        
            a = complex(i,j)
            a_abs = abs(a)
            
            if (a_abs <= alpha) and (a_abs not in mod_list):
            
                alpha_list.append(a)
                mod_list.append(a_abs)
    
    return alpha_list, mod_list


def Master_Equation(w0, beta_1, beta_2, tlist, L_operators, alpha):

    ## density matrices

    Z1 = 2*np.cosh(beta_1*w0/2)
    Z2 = 2*np.cosh(beta_2*w0/2)

    p1 = np.exp(beta_1*w0/2)/Z1
    p2 = np.exp(beta_2*w0/2)/Z2

    rho0_q1 = Qobj([[p1, 0],[0, 1-p1]])

    rho0_q2 = Qobj([[p2, 0],[0, 1-p2]])

    G = basis(2,0)
    E = basis(2,1)

    v00 = tensor(G, G)
    v01 = tensor(G, E)
    v10 = tensor(E, G)
    v11 = tensor(E, E)

    coherences_matrix = alpha * v01 * v10.dag() + alpha.conjugate() * v10 * v01.dag()

    rho0 = tensor(rho0_q1, rho0_q2) + coherences_matrix

    print(rho0.isherm, rho0.tr())
    
    evals, evecs = rho0.eigenstates()
    
    for val in evals:
        if val < 0:
            print(val)
    

    ## solve master equation

    dataME = mesolve(H_S, rho0, tlist, L_operators, [])

    rhof = dataME.states
    
    rhof_q1 = []
    rhof_q2 = []
    
    for rt in rhof:
    
        rhof_q1.append( rt.ptrace([0]) )
        rhof_q2.append( rt.ptrace([1]) )

    return rhof, rhof_q1, rhof_q2


def Write_Density_Matrices(rhof, rhof_q1, rhof_q2, path):

    ## write rho in output filess

    for t, rhot in enumerate(rhof):
        
        rhot_q1 = rhof_q1[t]
        rhot_q2 = rhof_q2[t]

        f = open(f'./{path}/rhof_t{t}.txt', 'w')
        f1 = open(f'./{path}/rhof_q1_t{t}.txt', 'w')
        f2 = open(f'./{path}/rhof_q2_t{t}.txt', 'w')

        for i in range(4):
            for j in range(4):
            
                f.write(f'{rhot.full()[i][j]} ')
                
            f.write('\n')
        
        for i in range(2):
            for j in range(2):
            
                f1.write(f'{rhot_q1.full()[i][j]} ')
                f2.write(f'{rhot_q2.full()[i][j]} ')

            f1.write('\n')
            f2.write('\n')
                
    f.close()
    f1.close()
    f2.close()

    


######################## MAIN ########################

## parameters

w0 = 1
gamma = 0.1
g = 0

tlist = np.arange(0, 20, 0.01)

Tf_banho = 1.0
Tf_qubit = 1.0

p_final = pFunc(Tf_qubit, w0)

Sr_inicial = 0.05


## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_int = tensor(sigmap(), sigmam()) + tensor(sigmam(), sigmap())

H_S = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2) + g*H_int


## temperaturas equidistantes

Tc_list = np.arange(0.05, 50, 0.0001)
Th_list = np.arange(0.05, 50, 0.0001)

pc, Src, Tc, Sr0, pList0 = Temperaturas_e_Populacoes(0, w0, Tc_list, p_final, Sr_inicial)

ph, Srh, Th, Sr1, pList1 = Temperaturas_e_Populacoes(1, w0, Th_list, p_final, Sr_inicial)

print(Tc, Tf_qubit, Th)


plt.plot(pList0, Sr0, color='orange', label=f'pf = {p_final:.3f}')
plt.plot(pList1, Sr1, color='orange')
plt.scatter([p_final], [0], color='orange', label=f'Tw = {Tf_qubit:.3f}')
plt.scatter([pc], [Src], color='blue', label=f'Tc = {Tc:.3f}')
plt.scatter([ph], [Srh], color='red', label=f'Th = {Th:.3f}')
plt.hlines(y=Sr_inicial, xmin=min(pList1), xmax=max(pList0), color='black', label='Initial Relative Entropy')
plt.xlabel('Populations')
plt.ylabel('Relative Entropy')
plt.legend()
plt.show()


beta_c = 1/Tc
beta_h = 1/Th
beta_R = 1/Tf_banho


## collapse operators

L_operators = []

evals, evecs = H_S.eigenstates()

for n, valn in enumerate(evals):
    for m, valm in enumerate(evals):
        
        dmn = valm.real - valn.real

        if dmn > 1e-9:
             
            f = Distribution('bose', beta_R, dmn)
                
            C_emi = np.sqrt(gamma * (1 + f)) * (evecs[n] * evecs[m].dag())
                
            L_operators.append(C_emi)
            
            if f > 0:
                C_abs = np.sqrt(gamma * f) * (evecs[m] * evecs[n].dag())

                L_operators.append(C_abs)


## master equation solver + write output file

os.mkdir(f'./DensityMatrices')


## Heating qubit 1 and Cooling qubit 2

f_alpha = open(f'./DensityMatrices/cmod.txt', 'w')

alpha_list, mod_list = Coherences(w0, beta_c, beta_h)

print(mod_list)

for a, alpha in enumerate(alpha_list):
    
    f_alpha.write(f'{mod_list[a]}\n')
    
    print('Heating 1 and Cooling 2', mod_list[a])

    ## master equation solver

    rhof, rhof_q1, rhof_q2 = Master_Equation(w0, beta_c, beta_h, tlist, L_operators, alpha)

    ## write output file
    
    path = f'./DensityMatrices/c_{mod_list[a]}'
    
    os.mkdir(path)
    
    Write_Density_Matrices(rhof, rhof_q1, rhof_q2, path)

f_alpha.close()
