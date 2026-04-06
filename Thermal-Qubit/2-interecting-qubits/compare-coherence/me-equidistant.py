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

    alpha = 1/(Z1*Z2)
    
    return alpha


def Collapse_Operators(beta_R, H_S, gamma):

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
    
    return L_operators



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

    dataME = mesolve(H_S, rho0, tlist, c_ops=L_operators, e_ops=[])

    rhof = dataME.states
    
    rhof_q1 = []
    rhof_q2 = []
    
    for rt in rhof:
    
        rhof_q1.append( rt.ptrace([0]) )
        rhof_q2.append( rt.ptrace([1]) )

    return rhof, rhof_q1, rhof_q2


def Write_Density_Matrices(rhof, rhof_q1, rhof_q2, c, g):

    ## write rho in output filess

    file_rhof = open(f'./DensityMatrices/rhof_c{c}_g{g}.txt', 'w')
    file_rhof_q1 = open(f'./DensityMatrices/rhof_q1_c{c}_g{g}.txt', 'w')
    file_rhof_q2 = open(f'./DensityMatrices/rhof_q2_c{c}_g{g}.txt', 'w')

    for t in range(len(rhof)):
        
        for i in range(4):
            for j in range(4):
        
                file_rhof.write(f'{rhof[t].full()[i][j]} ')
                
        file_rhof.write('\n')
        
        for i in range(2):
            for j in range(2):
        
                file_rhof_q1.write(f'{rhof_q1[t].full()[i][j]} ')
                file_rhof_q2.write(f'{rhof_q2[t].full()[i][j]} ')
                
        file_rhof_q1.write('\n')
        file_rhof_q2.write('\n')
        
    file_rhof.close()
    file_rhof_q1.close()
    file_rhof_q2.close()

    


######################## MAIN ########################

## parameters

w0 = 1
gamma = 0.1
g = 0.8

tlist = np.arange(0, 30, 0.01)

Tf_banho = 1.0
beta_R = 1/Tf_banho

Sr_inicial = 0.05


## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_int = tensor(sigmap(), sigmam()) + tensor(sigmam(), sigmap())

H_S = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2) + g*H_int


## qubits final temperature

L_operators = Collapse_Operators(beta_R, H_S, gamma)

rho_ss = steadystate(H_S, L_operators)

qubit = rho_ss.ptrace(0)

p0 = qubit.full()[0][0].real
p1 = qubit.full()[1][1].real

Tf_qubit = (w0/np.log(p0/p1))

p_final = pFunc(Tf_qubit, w0)

## temperaturas equidistantes

Tc_list = np.arange(0.05, 50, 0.0001)
Th_list = np.arange(0.05, 50, 0.0001)

pc, Src, Tc, Sr0, pList0 = Temperaturas_e_Populacoes(0, w0, Tc_list, p_final, Sr_inicial)

ph, Srh, Th, Sr1, pList1 = Temperaturas_e_Populacoes(1, w0, Th_list, p_final, Sr_inicial)

print(Tc, Tf_qubit, Th)


plt.scatter([p_final], [0], color='orange', label=f'Tw = {Tf_qubit:.3f}')
plt.scatter([pc], [Src], color='blue', label=f'Tc = {Tc:.3f}')
plt.scatter([ph], [Srh], color='red', label=f'Th = {Th:.3f}')
plt.plot(pList0, Sr0, color='orange')
plt.plot(pList1, Sr1, color='orange')
plt.hlines(y=Sr_inicial, xmin=min(pList1), xmax=max(pList0), color='black', label='Initial Relative Entropy')
plt.xlabel('Populations', fontsize=12)
plt.ylabel('Relative Entropy', fontsize=12)
plt.title(f'g = {g}', fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.legend(fontsize=12)
plt.show()


beta_c = 1/Tc
beta_h = 1/Th

fT = open('./DensityMatrices/temperature_qubit.txt', 'a')
fT.write(f'{g} {Tf_qubit} {Th} {Tc}\n')
fT.close()

## master equation solver + write output file


## Heating qubit 1 and Cooling qubit 2

alpha = Coherences(w0, beta_c, beta_h)

alpha_min_max = [0, alpha]
alphaName = ['min', 'max']

falpha = open('./DensityMatrices/coherences.txt', 'a')
falpha.write(f'{g} {alpha_min_max[0]} {alpha_min_max[1]}\n')
falpha.close()

for n, alpha in enumerate(alpha_min_max):
    
    print('Heating 1 and Cooling 2', alpha)

    ## master equation solver

    rhof, rhof_q1, rhof_q2 = Master_Equation(w0, beta_c, beta_h, tlist, L_operators, alpha)

    ## write output file
    
    Write_Density_Matrices(rhof, rhof_q1, rhof_q2, alphaName[n], g)
    

    
    
    
    
