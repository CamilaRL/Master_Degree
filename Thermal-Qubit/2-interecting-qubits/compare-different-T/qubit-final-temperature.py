import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
    
    
def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f
  
  
## parameters

w0 = 1
gamma = 0.1
glist = np.arange(0, 2, 0.1)

T_global_final = 1

beta_R = 1/T_global_final

dtol = 1e-9

## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H_int = tensor(sigmap(), sigmam()) + tensor(sigmam(), sigmap())


## lindabladian diagonalization

Tlist = []

for g in glist:
    ## system hamiltonian

    H_S = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2) + g*H_int

    ## collapse operators

    L_operators = []

    evals, evecs = H_S.eigenstates()
    
    for n, valn in enumerate(evals):
        for m, valm in enumerate(evals):
            
            dmn = valm.real - valn.real

            if dmn > dtol:
                 
                f = Distribution('bose', beta_R, dmn)
                    
                C_emi = np.sqrt(gamma * (1 + f)) * (evecs[n] * evecs[m].dag())

                L_operators.append(C_emi)
                
                if f > 0:
                    C_abs = np.sqrt(gamma * f) * (evecs[m] * evecs[n].dag())

                    L_operators.append(C_abs)
    
    
    rho_ss = steadystate(H_S, L_operators)

    qubit = rho_ss.ptrace(0)

    p0 = qubit.full()[0][0].real
    p1 = qubit.full()[1][1].real

    Tq = (w0/np.log(p0/p1))
    
    print(g, Tq)

    Tlist.append(Tq)
    

plt.plot(glist, Tlist)
plt.xlabel('g')
plt.ylabel('Temperature')
plt.title('Qubits final temperature')
plt.show()


