import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
import math
    
def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f


def Plot_Same_Delta(glist, mn_delta, deltalist, m, n):

    for i in range(len(deltalist)):
        plt.plot(glist, mn_delta[i], label=r'$\Delta$ ='+f'{deltalist[i]}')
    
    plt.xlabel('g')
    plt.ylabel('Eigenvalues Difference')
    plt.title(f'Transition from n = {n} to m = {m}')
    plt.legend()
    plt.show()
    

## parameters

w0 = 1
gamma = 0.1

T_global_final = 1

beta_R = 1/T_global_final

glist = np.arange(0, 1.1, 0.1)
J = 1
deltalist = [-1, 0, 0.1]

'''
## Coefficients in function of the eigenenergies difference

dmnList = np.arange(0.1, 5, 0.1)
fList = [Distribution('bose', beta_R, dmn) for dmn in dmnList]

c_emi = []
c_abs = []

for f in fList:

    c_emi.append(np.sqrt(gamma * (1 + f)))
    c_abs.append(np.sqrt(gamma * f))

plt.plot(dmnList, c_emi, label='Emission Coefficient')
plt.plot(dmnList, c_abs, label='Absorption Coefficient')
plt.xlabel('Eigenenergy Difference')
plt.ylabel('Coefficients')
plt.legend()
plt.show()
'''

## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H0 = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)

sxsx = tensor(sigmax(), sigmax())
sysy = tensor(sigmay(), sigmay())
szsz = tensor(sigmaz(), sigmaz())

## lindabladian diagonalization

m1_n0_delta = []
m2_n0_delta = []
m3_n0_delta = []
m2_n1_delta = []
m3_n1_delta = []
m3_n2_delta = []

for delta in deltalist:
    
    m1_n0 = []
    m2_n0 = []
    m3_n0 = []
    m2_n1 = []
    m3_n1 = []
    m3_n2 = []
    
    for g in glist:
        
        ## system hamiltonian

        H_int = (1 + g)*sxsx + (1 - g)*sysy + delta*szsz
        
        H_S = H0 + J*H_int
        
        evals, evecs = H_S.eigenstates()
        
        for n, valn in enumerate(evals):
            for m, valm in enumerate(evals):
                
                dmn = valm.real - valn.real
  
                if dmn > 1e-20:
                    
                    if m == 1 and n == 0:
                        m1_n0.append(dmn)
                        
                    elif m == 2 and n == 0:
                        m2_n0.append(dmn)
                    
                    elif m == 3 and n == 0:
                        m3_n0.append(dmn)
                    
                    elif m == 2 and n == 1:
                        m2_n1.append(dmn)
                    
                    elif m == 3 and n == 1:
                        m3_n1.append(dmn)
                    
                    elif m == 3 and n == 2:
                        m3_n2.append(dmn)
                    

    m1_n0_delta.append(m1_n0)
    m2_n0_delta.append(m2_n0)
    m3_n0_delta.append(m3_n0)
    m2_n1_delta.append(m2_n1)
    m3_n1_delta.append(m3_n1)
    m3_n2_delta.append(m3_n2)
                    

Plot_Same_Delta(glist, m1_n0_delta, deltalist, 1, 0)
Plot_Same_Delta(glist, m2_n0_delta, deltalist, 2, 0)
Plot_Same_Delta(glist, m3_n0_delta, deltalist, 3, 0)
Plot_Same_Delta(glist, m2_n1_delta, deltalist, 2, 1)
Plot_Same_Delta(glist, m3_n1_delta, deltalist, 3, 1)
Plot_Same_Delta(glist, m3_n2_delta, deltalist, 3, 2)



