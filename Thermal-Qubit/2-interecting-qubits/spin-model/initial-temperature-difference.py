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
            
    return pi


def Temperaturas_e_Populacoes(metade, w0, Tlist, p_final, Sr_inicial):
    
    pList = [pFunc(T, w0) for T in Tlist]
    
    p_i = Interseccao_Inicial(metade, p_final, pList, Sr_inicial)
    
    T_i = w0 / np.log(p_i/(1-p_i))

    return T_i

  
## parameters

w0 = 1
gamma = 0.1

Sr_inicial = 0.01
Tlist = np.arange(0.05, 25, 0.0001)

T_global_final = 1

beta_R = 1/T_global_final

dtol = 1e-9

glist = np.arange(-1, 1, 0.1)
J = 1
deltalist = [-1, -0.5, 0, 0.1]

## hamiltonians

H_q1 = -w0*sigmaz()/2

H_q2 = -w0*sigmaz()/2

H0 = tensor(H_q1, qeye(2)) + tensor(qeye(2), H_q2)

sxsx = tensor(sigmax(), sigmax())
sysy = tensor(sigmay(), sigmay())
szsz = tensor(sigmaz(), sigmaz())

## lindabladian diagonalization

for delta in deltalist:

    dT_list = []
    
    for g in glist:
    
        print(g)
        
        ## system hamiltonian

        H_int = (1 + g)*sxsx + (1 - g)*sysy + delta*szsz
        
        H_S = H0 + J*H_int 

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
        
        p_final = pFunc(Tq, w0)

        Tc = Temperaturas_e_Populacoes(0, w0, Tlist, p_final, Sr_inicial)

        Th = Temperaturas_e_Populacoes(1, w0, Tlist, p_final, Sr_inicial)

        dT_list.append(Th - Tc)
    
    plt.plot(glist, dT_list, linewidth=2, label=r'$\delta$ = '+f'{delta}')
    

plt.xlabel('g', fontsize=12)
plt.ylabel(r'$\Delta$T', fontsize=12)
plt.title('Qubits Initial Temperature Difference', fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.legend(loc='upper right', fontsize=12)
plt.show()


