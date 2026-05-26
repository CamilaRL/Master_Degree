import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import os
import math


def EoF(rho):

    rho = Qobj(rho, dims=[[2, 2], [2, 2]])
        
    C = concurrence(rho)
    
    C = max(0.0, min(1.0, C))
        
    x = 0.5*(1 + np.sqrt(1 - C**2))

    if x == 0 or x == 1:
        E = 0
            
    else:
        E = -x * math.log2(x) - (1-x) * math.log2(1-x)
        
    return E
    

def Rho_SS(beta_eq, H_S):
    
    rho_unnormalized = (-beta_eq * H_S).expm()
    
    Z = rho_unnormalized.tr()
    
    rho_ss = rho_unnormalized/Z
    
    return rho_ss
    




### MAIN ###

w0 = 1

beta_eq = 1

model_name = 'XXZ'
J = [-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1]
delta = [-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1]
gaminha = 0

size = len(J)

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_0 = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)

EoF_map = np.zeros((size, size))

for j in range(size):
    # eixo y
    for d in range(size):
        # eixo x
    
        H_int = (1+gaminha)*tensor(sigmax(), sigmax()) + (1-gaminha)*tensor(sigmay(), sigmay()) + delta[d]*tensor(sigmaz(), sigmaz())

        H_S = H_0 - J[j]*H_int
        
        rho_ss = Rho_SS(beta_eq, H_S)
        
        qubit = rho_ss.ptrace(0)

        p0 = qubit.full()[0][0].real
        p1 = qubit.full()[1][1].real

        Tf_qubit = (w0/np.log(p0/p1))

        p_final = np.exp(w0/(2*Tf_qubit))/(2*np.cosh(w0/(2*Tf_qubit)))
        
        EoF_map[j][d] = EoF(rho_ss)
        if Tf_qubit < 2 and EoF_map[j][d] > 0:
            print(J[j], delta[d], EoF_map[j][d], Tf_qubit)
    
    
image = plt.imshow(EoF_map, cmap='viridis', aspect='equal', origin='lower', vmin=0)

plt.yticks(ticks=np.arange(0, size, 1), labels=J)
plt.xticks(ticks=np.arange(0, size, 1), labels=delta)
plt.colorbar()
plt.title(f'Entanglement of Formation | {model_name} Model')
plt.ylabel('J')
plt.xlabel(r'$\Delta$')
plt.tight_layout()
plt.show()