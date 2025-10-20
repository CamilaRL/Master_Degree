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

def Write_Outfile(temperature, fisher, path):

    f = open(path, 'w')
    
    for i in range(len(temperature)):
        
        f.write(f'{temperature[i]} {fisher[i]}\n')
        
    f.close()
    

# informacao de fisher sistema inteiro
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

dtT = tT - 1
ddTemp = linspace(TempMin,TempMax,dtT) # Temperature

dTemp = Temp[1] - Temp[0]


## time

tSEmax = 10*J
tSE = np.arange(0.001, tSEmax, 0.001) # Time S-E


dtSE = tSE[1]-tSE[0]


n = 30
NN = range(0, n)
Tp = len(Temp)

###############################################################################

# POVM Parameter
theta = linspace(0,pi,10)
phi = linspace(0,2*pi,20)

dTheta = theta[1]-theta[0]
dPhi = phi[1]-phi[0]

#################


nthermo = np.zeros((len(Temp)))

FthGab = np.zeros((len(Temp)))


p1dA = np.zeros((len(Temp),len(theta),len(phi)))
p2dA = np.zeros((len(Temp),len(theta),len(phi)))

p1lndA = np.zeros((len(Temp),len(theta),len(phi)))
p2lndA = np.zeros((len(Temp),len(theta),len(phi)))

derP1A = np.zeros((len(Temp)-1,len(theta),len(phi)))
derP2A = np.zeros((len(Temp)-1,len(theta),len(phi)))


#################
F1A =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA = np.zeros((len(Temp)-1,len(theta)))


QFIA = np.zeros((len(Temp)-1,len(range(0, n))))


QFI_Therm =  np.zeros((len(Temp)))


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
HP4B = -g * (Sx1 + Sx2 + Sx3 + Sx4 + Sx5)

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

S01 = tensor(A1,qeye(2),qeye(2),qeye(2),qeye(2))
S02 = tensor(qeye(2),A1,qeye(2),qeye(2),qeye(2))
S03 = tensor(qeye(2),qeye(2),A1,qeye(2),qeye(2))
S04 = tensor(qeye(2),qeye(2),qeye(2),A1,qeye(2))
S05 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),A1)

S0 = S01 + S02 + S03 + S04 + S05

Sm = Sm1 + Sm2 + Sm3 + Sm4 + Sm5
Sp = Sp1 + Sp2 + Sp3 + Sp4 + Sp5

###############################################################################
r = -1

for T in Temp:

    r = r + 1
    print(r)
###############################################################################
###############################  Thermal Number  ##############################

    nt = 1/(exp((OmegaS)/(T))-1)
    nthermo[r] = nt


    FthGab[r] = ((OmegaS/(T**2))**2)*(1/(np.cosh(OmegaS/T)))**2
###############################################################################
###############################  Collapse Operator  ###########################

    C1S = np.sqrt(gamma*(nt+1))*Sm
    C2S = np.sqrt(gamma*nt)*Sp

    Clist = [C1S,C2S] 

###############################################################################
#################################  Initial State  #############################

    G = basis(2,0) # base: excited state
    E = basis(2,1) # base: ground state


#######  Ancilla 1
    A1_st = G

    A1_mat = A1_st*A1_st.dag()

#######  Ancilla 2
    A2_st = G

    A2_mat = A2_st*A2_st.dag()

#######  Ancilla 3
    A3_st = G

    A3_mat = A3_st*A3_st.dag()

 #######  Ancilla 4   
    A4_st = G

    A4_mat = A4_st*A4_st.dag()

 #######  Ancilla 5
    A5_st = G

    A5_mat = A5_st*A5_st.dag()

    
    psi = tensor(A1_st,A2_st,A3_st,A4_st,A5_st)

#######  Completo

    S0mat = psi*psi.dag()

    medataSE = mesolve(H,S0mat,tSE,[Clist],[]) # Master equation evolution - Sistema + Ambiente
    rho0 = medataSE.states # Take matrices in each time
    medataH1 = mesolve(H,S0mat,tSE,[Clist],[Sz1]) # Master equation evolution - Sistema + Ambiente
    Exp_rho0 = medataH1.expect[0] # Take matrices in each time

    
    rhof = rho0[-1] # Take matrices in each time
    Rhof1 = rho0[-1].ptrace([0])
    Rhof2 = rho0[-1].ptrace([1])
    Rhof3 = rho0[-1].ptrace([2])
    Rhof4 = rho0[-1].ptrace([3])
    Rhof5 = rho0[-1].ptrace([4])

    
    medataH2 = mesolve(H,S0mat,tSE,[Clist],[Sz1*Sz1]) # Master equation evolution - Sistema + Ambiente
    Exp_rho02 = medataH2.expect[0]

    

    Therm1 = (Exp_rho02[-1] - (Exp_rho0[-1]*Exp_rho0[-1]))


    QFI_Therm[r] = Therm1/(T**4)
    

        
###############################################################################    
###############################   POVM   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1

            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()
            P1 = tensor(P1a,P1a,P1a,P1a,P1a)
            

            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()
            P2 = tensor(P2a,P2a,P2a,P2a,P2a)


            p1A = (P1*rhof).tr()

            R1A = P1*rhof

            p2A = (P2*rhof).tr()
            R2A = P2*rhof
        
            p1dA[r,q,z] = real(p1A)
            p2dA[r,q,z] = real(p2A)
        
            p1lnA = np.log(p1A)
            p2lnA = np.log(p2A)
            
            p1lndA[r,q,z] = real(p1lnA)
            p2lndA[r,q,z] = real(p2lnA)

                
###############################################################################
##################################   Fisher Total  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):
                # print(i)
        for p in range(len(phi)):
            derP1A[r,i,p] = (p1lndA[r+1,i,p] - p1lndA[r,i,p])/dtSE
            derP2A[r,i,p] = (p2lndA[r+1,i,p] - p2lndA[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 
    for i in range(len(theta)):
                # print(i)
        for p in range(len(phi)):
            f1 = p1dA[r,i,p]*(derP1A[r,i,p])**2
            f2 = p2dA[r,i,p]*(derP2A[r,i,p])**2
                    # print(f1)
            fx = f1+f2
                    # FI = max(fx)
    
            F1A[r,i,p] = f1
            F2A[r,i,p] = f2
    
            FxA[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA[r,i] = max(FxA[r,i,:])

for r in range(len(Temp)-1):
        QFIA[r] = max(QFI_tA[r,:])
            

###############################################################################
FthSca = np.zeros((len(Temp)-1))
for r in range(len(Temp)-1):
    
    FthSca[r] = (1/(nthermo[r+1]*(nthermo[r+1]+1)*(2*nthermo[r+1]+1)**2))*((nthermo[r+1]-nthermo[r])/dTemp)**2

###############################################################################
Write_Outfile(ddTemp, QFIA[:,-1], f'./Results/QFI_all_g{g:.2f}_ttherm{tSEmax:.3f}.txt')



tend = time.time()    # tempo final de processamento
delta = tend - tin    # funcao calculo do intervalo de tempo  de processamento
print (delta)


