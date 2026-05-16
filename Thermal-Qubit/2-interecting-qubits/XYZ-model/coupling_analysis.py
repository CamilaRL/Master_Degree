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
            Ti = w0 / np.log(pi/(1-pi))
            
    return pi, Ti, Sri
    

def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f


def Coherences(w0, beta_1, beta_2, rho0):

    alpha_list = []
    mod_list = []

    Z1 = 2*np.cosh(beta_1*w0/2)
    Z2 = 2*np.cosh(beta_2*w0/2)

    alpha = 1/(Z1*Z2)
    
    G = basis(2,0)
    E = basis(2,1)

    v00 = tensor(G, G)
    v01 = tensor(G, E)
    v10 = tensor(E, G)
    v11 = tensor(E, E)
    
    for f in np.arange(0, alpha, 0.001):
        for e in np.arange(0, f+0.001, 0.001):
            
            coherences_matrix = e * v01 * v10.dag() + e.conjugate() * v10 * v01.dag()
            entanglement_matrix = f * v00 * v11.dag() + f.conjugate() * v11 * v00.dag()

            rho = rho0 + coherences_matrix + entanglement_matrix
            
            C = concurrence(Qobj(rho, dims=[[2, 2], [2, 2]]))
            
            sysy = tensor(sigmay(), sigmay())
            
            rho_rhotilde = rho * sysy * rho * sysy
            
            evals = rho_rhotilde.eigenenergies()
            
            evals = sorted(evals)

            evals_diff = np.sqrt(evals[3]) - np.sqrt(evals[2]) - np.sqrt(evals[1]) - np.sqrt(evals[0])
            
            
            lamb1 = alpha + e
            lamb2 = alpha - e
            lamb3 = alpha + f
            lamb4 = alpha - f
                
            lambList = sorted([lamb1, lamb2, lamb3, lamb4])

            lamb_diff = lambList[3] - lambList[2] - lambList[1] - lambList[0]
            
            if lamb_diff > 0:
                print(alpha, f, e, lamb_diff)
                #print(C, max(0, evals_diff), max(0, lamb_diff))
    
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


### MAIN ###


## parameters

w0 = 1
gamma = 0.1
JList = np.arange(0, 0.8, 0.1)
gaminha = 0
delta = 1

Tf_banho = 1.0
beta_R = 1/Tf_banho

Sr_inicial = 0.05

TList = np.arange(0.05, 50, 0.0001)
pList = [pFunc(T, w0) for T in TList]


## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_0 = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)

H_int = (1+gaminha)*tensor(sigmax(), sigmax()) + (1-gaminha)*tensor(sigmay(), sigmay()) + delta*tensor(sigmaz(), sigmaz())


Tf_list = []
Tc_list = []
Th_list = []
alpha_list = []

for J in JList:
    print(J)
    H_S = H_0 - J*H_int

    ## qubits final temperature

    L_operators = Collapse_Operators(beta_R, H_S, gamma)

    rho_ss = steadystate(H_S, L_operators)

    qubit = rho_ss.ptrace(0)

    p0 = qubit.full()[0][0].real
    p1 = qubit.full()[1][1].real

    Tf_qubit = (w0/np.log(p0/p1))

    p_final = pFunc(Tf_qubit, w0)


    ## temperaturas equidistantes

    pc, Tc, Src = Interseccao_Inicial(0, p_final, pList, Sr_inicial)
    ph, Th, Srh = Interseccao_Inicial(1, p_final, pList, Sr_inicial)

    beta_c = 1/Tc
    beta_h = 1/Th
    
    Tf_list.append(Tf_qubit)
    Tc_list.append(Tc)
    Th_list.append(Th)
    
    Z1 = 2*np.cosh(beta_c*w0/2)
    Z2 = 2*np.cosh(beta_h*w0/2)

    p1 = np.exp(beta_c*w0/2)/Z1
    p2 = np.exp(beta_h*w0/2)/Z2

    rho0_q1 = Qobj([[p1, 0],[0, 1-p1]])

    rho0_q2 = Qobj([[p2, 0],[0, 1-p2]])
    
    rho0 = tensor(rho0_q1, rho0_q2)
    
    alpha = Coherences(w0, beta_c, beta_h, rho0)
    
    alpha_list.append(alpha)



plt.plot(JList, Tf_list, color='orange', label='Final Temperature')
plt.plot(JList, Tc_list, color='Blue', label='Colder Initial Temperature')
plt.plot(JList, Th_list, color='Red', label='Hotter Initial Temperature')
plt.ylabel('Temperature')
plt.xlabel('J')
plt.legend()
plt.show()


plt.plot(JList, alpha_list, color='black')
plt.ylabel(r'$\alpha$')
plt.xlabel('J')
plt.show()