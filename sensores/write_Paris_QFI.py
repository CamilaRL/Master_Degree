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
    

def Quantum_Fisher_Information(rho_list, temp_list):

    
    ### derivada de rho em relacao a temperatura
    
    size_rhoT = rho_list[0].shape[0]
    
    drho_list = [Qobj(np.zeros((size_rhoT, size_rhoT))) for T in range(len(temp_list)-1)]
    
    for i in range(size_rhoT):
        for j in range(size_rhoT):
    
            for t in range(len(temp_list)-1):
    
                rho_i = rho_list[t][i][j]
                rho_f = rho_list[t+1][i][j]
                temp_i = temp_list[t]
                temp_f = temp_list[t+1]
                
                drho_list[t][i][j] = (rho_f - rho_i)/(temp_f - temp_i)
    
    
    ### calculo da Infomacao de Fisher Quantica
    QFI_T = []
    
    for t in range(len(temp_list)-1):
    
        autoval, autovec = rho_list[t].eigenstates()
        
        QFI = 0
        
        print(autovec[0].dag() * drho_list[t])
        '''for n in range(len(autoval)):
            for m in range(len(autoval)):
            
                #QFI = QFI + (abs(autovec[n].dag() * drho_list[t] * autovec[m])**2)/(autoval[n] + autoval[m])
                print(autovec[n].dag() * drho_list[t] * autovec[m])
        print(QFI)
        QFI_T.append(QFI)'''
        
    

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

g = 1.0  # Sistema-Campo
J = 1.0       # Sistema-sistema

## temperature

tT = 100
TempMin = 0.005
TempMax = 2.0
Temp = linspace(TempMin,TempMax,tT) # Temperature

dtT = tT - 1
ddTemp = linspace(TempMin,TempMax,dtT) # Temperature

dTemp = Temp[1] - Temp[0]


## time

tp = 100 # Step
tSEmax = (pi/2) -1
tSE = linspace(0.00001,tSEmax,tp)  # Time S-E

td = linspace(0.0,5.0,tp-1)

dtSE = tSE[1]-tSE[0]


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

# System: 4 qubits
# qubit 1:
Sm1 = tensor(sigmap(),qeye(2),qeye(2),qeye(2))
Sp1 = tensor(sigmam(),qeye(2),qeye(2),qeye(2))
Sx1 = tensor(sigmax(),qeye(2),qeye(2),qeye(2))
Sy1 = tensor(sigmay(),qeye(2),qeye(2),qeye(2))
Sz1 = tensor(sigmaz(),qeye(2),qeye(2),qeye(2))

# print(Sp.ptrace(0))

# qubit 2:
Sm2 = tensor(qeye(2),sigmap(),qeye(2),qeye(2))
Sp2 = tensor(qeye(2),sigmam(),qeye(2),qeye(2))
Sx2 = tensor(qeye(2),sigmax(),qeye(2),qeye(2))
Sy2 = tensor(qeye(2),sigmay(),qeye(2),qeye(2))
Sz2 = tensor(qeye(2),sigmaz(),qeye(2),qeye(2))

# qubit 3:
Sm3 = tensor(qeye(2),qeye(2),sigmap(),qeye(2))
Sp3 = tensor(qeye(2),qeye(2),sigmam(),qeye(2))
Sx3 = tensor(qeye(2),qeye(2),sigmax(),qeye(2))
Sy3 = tensor(qeye(2),qeye(2),sigmay(),qeye(2))
Sz3 = tensor(qeye(2),qeye(2),sigmaz(),qeye(2))

# qubit 4:
Sm4 = tensor(qeye(2),qeye(2),qeye(2),sigmap())
Sp4 = tensor(qeye(2),qeye(2),qeye(2),sigmam())
Sx4 = tensor(qeye(2),qeye(2),qeye(2),sigmax())
Sy4 = tensor(qeye(2),qeye(2),qeye(2),sigmay())
Sz4 = tensor(qeye(2),qeye(2),qeye(2),sigmaz())


###############################################################################
###################################  Hamiltonian ##############################

HP4A = -J * (Sz1* Sz2 + Sz2* Sz3+ Sz3* Sz4)
HP4B = -g * (Sx1 + Sx2 + Sx3 + Sx4)

HK4A = -J * (Sz1* (Sz2 + Sz3 + Sz4) + Sz2* (Sz3 +Sz4) + Sz3* Sz4)
HK4B = -g * (Sx1 + Sx2 + Sx3 + Sx4)

HSd4A = -J * (Sz1* (Sz2 + Sz3 + Sz4) + Sz2* (Sz3) + Sz3* Sz4)
HSd4B = -g * (Sx1 + Sx2 + Sx3 + Sx4)

HC4A = -J * (Sz1* Sz2 + Sz2* Sz3+ Sz3* Sz4 + Sz4*Sz1)
HC4B = -g * (Sx1 + Sx2 + Sx3 + Sx4)

HpanA = -J * (Sz1* Sz2 + Sz2* Sz3+ Sz3* Sz4 + Sz4*Sz2)
HpanB = -g * (Sx1 + Sx2 + Sx3 + Sx4)

HS3A = -J * (Sz1* (Sz2 + Sz3 + Sz4))
HS3B = -g * (Sx1 + Sx2 + Sx3 + Sx4)

H = HP4A + HP4B
# H = HC4A + HC4B
# H = HK4A + HK4B
# H = HSd4A + HSd4B
# H = HpanA + HpanB
# H = HS3A + HS3B


G = basis(2,0) # base: excited state
E = basis(2,1) # base: ground state

A1 = G*G.dag()
A2 = E*E.dag()

S01 = tensor(A1,qeye(2),qeye(2),qeye(2))
S02 = tensor(qeye(2),A1,qeye(2),qeye(2))
S03 = tensor(qeye(2),qeye(2),A1,qeye(2))
S04 = tensor(qeye(2),qeye(2),qeye(2),A1)

S01 = tensor(A1,qeye(2),qeye(2),qeye(2))
S02 = tensor(qeye(2),A1,qeye(2),qeye(2))
S03 = tensor(qeye(2),qeye(2),A1,qeye(2))
S04 = tensor(qeye(2),qeye(2),qeye(2),A1)

S0 = S01 + S02 + S03 + S04

Sm = Sm1 + Sm2 + Sm3 + Sm4
Sp = Sp1 + Sp2 + Sp3 + Sp4

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

#######  Completo

    psi = tensor(A1_st,A2_st,A3_st,A4_st)
    
    S0mat = psi*psi.dag()

#######  Dinamica

    medataSE = mesolve(H,S0mat,tSE,[Clist],[]) # Master equation evolution - Sistema + Ambiente
    rho0 = medataSE.states # Take matrices in each time
    
    rhof.append(rho0[-1]) # Take matrices in each time
    rhof1.append(rho0[-1].ptrace([0]))
    rhof2.append(rho0[-1].ptrace([1]))
    rhof3.append(rho0[-1].ptrace([2]))
    rhof4.append(rho0[-1].ptrace([3]))


################################# Quantum Fisher Information ########################################

Quantum_Fisher_Information(rhof, Temp)

            
###############################################################################
#Write_Outfile(ddTemp, QFIA1[:,-1], f'./Results/QFI_all_g{g:.1f}_ttherm{tSEmax:.3f}.txt')



tend = time.time()    # tempo final de processamento
delta = tend - tin    # funcao calculo do intervalo de tempo  de processamento
print (delta)


