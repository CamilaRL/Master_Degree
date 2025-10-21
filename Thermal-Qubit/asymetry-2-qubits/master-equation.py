import numpy as np
import matplotlib.pyplot as plt
from qutip import *


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


## parameters

w0 = 1
Ts = 2
Tr = 10

gamma = 1

beta_1 = 1/Ts
beta_2 = 1/Ts
beta_R = 1/Tr 

tlist = np.arange(0.005, 5, 0.001)

## hamiltonians

H_q1 = w0*sigmaz()/2

H_q2 = w0*sigmaz()/2

H_S = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)


## density matrices

Z1 = 2*np.cosh(beta_1*w0/Ts)
Z2 = 2*np.cosh(beta_2*w0/Ts)

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

coherences_matrix = alpha_list[0] * v01 * v10.dag() + alpha_list[0].conjugate() * v10 * v01.dag()

rho0 = tensor(rho0_q1, rho0_q2) + coherences_matrix


## collapse operators

L_operators = []

evals, evecs = H_S.eigenstates()

for n, valn in enumerate(evals):
    for m, valm in enumerate(evals):
        
        dmn = valm - valn
        
        if dmn < 0 : ## absortion
            
            f = Distribution('bose', beta_R, -dmn)
            
            C_abs = gamma * np.sqrt(1 + f) * (evecs[n] * evecs[m].dag())

            L_operators.append(C_abs)

        elif dmn > 0: ## emission
            
            f = Distribution('bose', beta_R, dmn)
                
            C_emi = gamma * np.sqrt(f) * (evecs[m] * evecs[n].dag())

            L_operators.append(C_emi)
            

## solve master equation

dataME = mesolve(H_S, rho0, tlist, L_operators, [])

rhof = dataME.states

## write rho in output filess

for t, rhot in enumerate(rhof):

    f = open(f'./DensityMatrices/rhof_t{t}.txt', 'w')

    for i in range(4):
        for j in range(4):
        
            f.write(f'{rhot[i][j]} ')
           
        f.write('\n')
            
f.close()







