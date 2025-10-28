from qutip import *
from math import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import time
from numpy import linalg as LA
from numpy.linalg import norm
import scipy.constants as constant

tin = time.time()

###############################################################################
#################################  FUNCTIONS  ################################
###############################################################################

def Write_Outfile(temperature, fisher, path):

    f = open(path, 'w')
    
    for i in range(len(temperature)):
        
        f.write(f'{temperature[i]} {fisher[i]}\n')
        
    f.close()
    

def Quantum_Fisher_Information(rho_list, temp_list, dimensao):

    
    ### derivada de rho em relacao a temperatura
    
    size_rhoT = rho_list[0].shape[0]
    
    drho_list = [np.zeros((size_rhoT, size_rhoT), dtype=complex) for k in range(len(temp_list)-1)]
    
    for i in range(size_rhoT):
        for j in range(size_rhoT):
            
            for t in range(len(temp_list)-1):
    
                rho_i = rho_list[t].full()[i][j]
                rho_f = rho_list[t+1].full()[i][j]
                temp_i = temp_list[t]
                temp_f = temp_list[t+1]

                drho_list[t][i][j] = (rho_f - rho_i)/(temp_f - temp_i)
    
    ### calculo da Infomacao de Fisher Quantica
    QFI_T = []
    
    for t in range(len(temp_list)-1):
        
        
        autoval, autovec = rho_list[t].eigenstates()
        
        QFI = complex(0,0)       
        
        for n in range(len(autoval)):
            for m in range(len(autoval)):
                
                bra = autovec[m].dag()
                ket = autovec[n]
                
                bra_drho_ket = (bra * Qobj(drho_list[t], dims=dimensao) * ket)

                QFI = QFI + (2*(bra_drho_ket.norm())**2)/(abs(autoval[n]) + abs(autoval[m]))
        
        QFI_T.append(QFI.real)
        
    return temp_list[:-1], QFI_T

###############################################################################
#################################  MAIN  ################################
###############################################################################




###############################################################################
#################################  Parameters  ################################
# hp = 6.62607015E-34 #m2 kg / s
# kb = 1.380649E-23   #m2 kg s-2 K-1

gamma = 1.0  # System decay rate
gammaE = 1.0  # Environment decay rate


OmegaS = 1.0#1E10*(T/2)  # System: energy level splitting
Omega1 = 1.0
Omega2 = 1.0
Omega3 = 1.0
Omega4 = 1.0
Omega5 = 1.0
Omega6 = 1.0

g = 0.45  # Sistema-Campo
J = 1.0       # Sistema-sistema

## temperature

tT = 100
TempMin = 0.005
TempMax = 2.0
Temp = linspace(TempMin,TempMax,tT) # Temperature


## time

tSEmax = 10*J
tSE = np.arange(0.001, tSEmax, 0.001) # Time S-E


n = 30
NN = range(0, n)
Tp = len(Temp)

###############################################################################
###################### Inicializacao Variaveis e Listas ######################
rhof = []
rhof1 = []
rhof2 = []
rhof3 = []
rhof4 = []
rhof5 = []


nthermo = np.zeros((len(Temp)))

FthGab = np.zeros((len(Temp)))



###############################################################################
#################################  Operadores  ################################

N=2

# System Alone:
sm= sigmap()
sp= sigmam()
sx= sigmax()
sy= sigmay()
sz= sigmaz()

# System: 5 qubits
# qubit 1:
Sm1 = tensor(sigmap(),qeye(2),qeye(2),qeye(2),qeye(2))
Sp1 = tensor(sigmam(),qeye(2),qeye(2),qeye(2),qeye(2))
Sx1 = tensor(sigmax(),qeye(2),qeye(2),qeye(2),qeye(2))
Sy1 = tensor(sigmay(),qeye(2),qeye(2),qeye(2),qeye(2))
Sz1 = tensor(sigmaz(),qeye(2),qeye(2),qeye(2),qeye(2))

# print(Sp.ptrace(0))

# qubit 2:
Sm2 = tensor(qeye(2),sigmap(),qeye(2),qeye(2),qeye(2))
Sp2 = tensor(qeye(2),sigmam(),qeye(2),qeye(2),qeye(2))
Sx2 = tensor(qeye(2),sigmax(),qeye(2),qeye(2),qeye(2))
Sy2 = tensor(qeye(2),sigmay(),qeye(2),qeye(2),qeye(2))
Sz2 = tensor(qeye(2),sigmaz(),qeye(2),qeye(2),qeye(2))

# qubit 3:
Sm3 = tensor(qeye(2),qeye(2),sigmap(),qeye(2),qeye(2))
Sp3 = tensor(qeye(2),qeye(2),sigmam(),qeye(2),qeye(2))
Sx3 = tensor(qeye(2),qeye(2),sigmax(),qeye(2),qeye(2))
Sy3 = tensor(qeye(2),qeye(2),sigmay(),qeye(2),qeye(2))
Sz3 = tensor(qeye(2),qeye(2),sigmaz(),qeye(2),qeye(2))

# qubit 4:
Sm4 = tensor(qeye(2),qeye(2),qeye(2),sigmap(),qeye(2))
Sp4 = tensor(qeye(2),qeye(2),qeye(2),sigmam(),qeye(2))
Sx4 = tensor(qeye(2),qeye(2),qeye(2),sigmax(),qeye(2))
Sy4 = tensor(qeye(2),qeye(2),qeye(2),sigmay(),qeye(2))
Sz4 = tensor(qeye(2),qeye(2),qeye(2),sigmaz(),qeye(2))

# qubit 5:
Sm5 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmap())
Sp5 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmam())
Sx5 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmax())
Sy5 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmay())
Sz5 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmaz())


###############################################################################
###################################  Hamiltonian ##############################

HP4A = -J * (Sz1* Sz2 + Sz2* Sz3 + Sz3* Sz4 + Sz4* Sz5)
HP4B = -g * (Sx1 + Sx2 + Sx3 + Sx4 + Sz5)

H = HP4A + HP4B

G = basis(2,0) # base: excited state
E = basis(2,1) # base: ground state

A1 = G*G.dag()
A2 = E*E.dag()

S01 = tensor(A1,qeye(2),qeye(2),qeye(2),qeye(2))
S02 = tensor(qeye(2),A1,qeye(2),qeye(2),qeye(2))
S03 = tensor(qeye(2),qeye(2),A1,qeye(2),qeye(2))
S04 = tensor(qeye(2),qeye(2),qeye(2),A1,qeye(2))
S05 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),A1)


S0 = S01 + S02 + S03 + S04 + S05

Sm = Sm1 + Sm2 + Sm3 + Sm4 + Sm5 
Sp = Sp1 + Sp2 + Sp3 + Sp4 + Sp5

###############################################################################

for r, T in enumerate(Temp):

    print(r)
###############################################################################
###############################  Thermal Number  ##############################

    nt = 1/(exp((OmegaS)/(T))-1)
    nthermo[r] = nt

###############################################################################
###############################  Collapse Operator  ###########################

    C1S = np.sqrt(gamma*(nt+1))*Sm
    C2S = np.sqrt(gamma*nt)*Sp

    Clist = [C1S,C2S] 

###############################################################################
#################################  Initial State  #############################

    G = basis(2,0) # base: excited state
    E = basis(2,1) # base: ground state


#######  qubit 1
    A1_st = G

    A1_mat = A1_st*A1_st.dag()

#######  qubit 2
    A2_st = G

    A2_mat = A2_st*A2_st.dag()

#######  qubit 3
    A3_st = G

    A3_mat = A3_st*A3_st.dag()

 #######  qubit 4   
    A4_st = G

    A4_mat = A4_st*A4_st.dag()
    
 #######  qubit 5
    A5_st = G

    A5_mat = A5_st*A5_st.dag()

#######  Completo

    psi = tensor(A1_st,A2_st,A3_st,A4_st,A5_st)
    
    S0mat = psi*psi.dag()

#######  Dinamica

    medataSE = mesolve(H,S0mat,tSE,[Clist],[]) # Master equation evolution - Sistema + Ambiente
    rho0 = medataSE.states # Take matrices in each time
    
    rhof.append(rho0[-1]) # Take matrices in each time
    rhof1.append(rho0[-1].ptrace([0]))
    rhof2.append(rho0[-1].ptrace([1]))
    rhof3.append(rho0[-1].ptrace([2]))
    rhof4.append(rho0[-1].ptrace([3]))
    rhof5.append(rho0[-1].ptrace([4]))

################################# Quantum Fisher Information ########################################

temp_all, QFI_all = Quantum_Fisher_Information(rhof, Temp, [[2,2,2,2,2],[2,2,2,2,2]])
temp_1, QFI_1 = Quantum_Fisher_Information(rhof1, Temp, [[2],[2]])
temp_2, QFI_2 = Quantum_Fisher_Information(rhof2, Temp, [[2],[2]])
temp_3, QFI_3 = Quantum_Fisher_Information(rhof3, Temp, [[2],[2]])
temp_4, QFI_4 = Quantum_Fisher_Information(rhof4, Temp, [[2],[2]])
temp_5, QFI_5 = Quantum_Fisher_Information(rhof5, Temp, [[2],[2]])


###############################################################################
Write_Outfile(temp_all, QFI_all, f'./Results/QFI_all_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(temp_1, QFI_1, f'./Results/QFI_1_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(temp_2, QFI_2, f'./Results/QFI_2_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(temp_3, QFI_3, f'./Results/QFI_3_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(temp_4, QFI_4, f'./Results/QFI_4_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(temp_5, QFI_5, f'./Results/QFI_5_g{g:.2f}_ttherm{tSEmax:.3f}.txt')


tend = time.time()    # tempo final de processamento
delta = tend - tin    # funcao calculo do intervalo de tempo  de processamento
print (delta)


