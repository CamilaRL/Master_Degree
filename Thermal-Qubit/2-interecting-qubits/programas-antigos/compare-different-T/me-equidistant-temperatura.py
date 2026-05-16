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
    

def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f


def Coherences(w0, beta_1, beta_2):

    alpha_list = []
    mod_list = []

    Z1 = 2*np.cosh(beta_1*w0)
    Z2 = 2*np.cosh(beta_2*w0)

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
g = 0.8

tlist = np.arange(0.1, 20, 0.01)

Tf_banho = 1.0
Tf_qubit = 1.1542306681364487

p_final = pFunc(Tf_qubit, w0)

## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_int = tensor(sigmap(), sigmam()) + tensor(sigmam(), sigmap())

H_S = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2) + g*H_int


## temperaturas equidistantes do caso sem acoplamento


Tc = 0.47144190368115213
Th = 5.668454511145365

beta_c = 1/Tc
beta_h = 1/Th
beta_R = 1/Tf_banho

pc = pFunc(Tc, w0)
ph = pFunc(Th, w0)

Sri_c = Entropia_Relativa_Populacoes(pc, p_final)
Sri_h = Entropia_Relativa_Populacoes(ph, p_final)

Tlist = np.arange(0.05, 50, 0.0001)
pList = [pFunc(T, w0) for T in Tlist]
Sr = [Entropia_Relativa_Populacoes(p, p_final) for p in pList]

plt.plot(pList, Sr, color='orange', label=f'pf = {p_final:.3f}')
plt.scatter([p_final], [0], color='orange', label=f'Tw = {Tf_qubit:.3f}')
plt.scatter([pc], [Sri_c], color='blue', label=f'Tc = {Tc:.3f}')
plt.scatter([ph], [Sri_h], color='red', label=f'Th = {Th:.3f}')
plt.hlines(y=Sri_c, xmin=min(pList), xmax=max(pList), color='blue', linestyle='--')
plt.hlines(y=Sri_h, xmin=min(pList), xmax=max(pList), color='red', linestyle='--')
plt.xlabel('Populations')
plt.ylabel('Relative Entropy')
plt.legend()
plt.show()


'''
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

for a, alpha in enumerate(alpha_list):
    
    f_alpha.write(f'{mod_list[a]}\n')
    
    print('Heating 1 and Cooling 2', a)

    ## master equation solver

    rhof, rhof_q1, rhof_q2 = Master_Equation(w0, beta_c, beta_h, tlist, L_operators, alpha)

    ## write output file
    
    path = f'./DensityMatrices/c_{mod_list[a]}'
    
    os.mkdir(path)
    
    Write_Density_Matrices(rhof, rhof_q1, rhof_q2, path)

f_alpha.close()





'''