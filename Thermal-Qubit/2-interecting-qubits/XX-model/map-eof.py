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

model_name = 'XX'
J = np.arange(-1, 1.1, 0.1)
delta = 0
gaminha = 0

size = len(J)

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_0 = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)

EoF_map = []

for j in range(size):
    # eixo x
    
    H_int = (1+gaminha)*tensor(sigmax(), sigmax()) + (1-gaminha)*tensor(sigmay(), sigmay()) + delta*tensor(sigmaz(), sigmaz())

    H_S = H_0 - J[j]*H_int
        
    rho_ss = Rho_SS(beta_eq, H_S)
            
    qubit = rho_ss.ptrace(0)

    p0 = qubit.full()[0][0].real
    p1 = qubit.full()[1][1].real

    Tf_qubit = (w0/np.log(p0/p1))

    p_final = np.exp(w0/(2*Tf_qubit))/(2*np.cosh(w0/(2*Tf_qubit)))
            
    EoF_map.append(EoF(rho_ss))
        
    if Tf_qubit < 2 and EoF_map[j] > 0:
        print(J[j], EoF_map[j], Tf_qubit)
    
    
plt.plot(J, EoF_map, color='black')
plt.title(f'Entanglement of Formation | {model_name} Model')
plt.ylabel('EoF')
plt.xlabel('J')
plt.show()