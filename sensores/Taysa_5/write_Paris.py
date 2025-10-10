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



###############################################################################
#################################  Parameters  ################################

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


ES = np.zeros((len(Temp),len(range(0, n))))


nthermo = np.zeros((len(Temp)))

FthGab = np.zeros((len(Temp)))


p1dA1 = np.zeros((len(Temp),len(theta),len(phi)))
p2dA1 = np.zeros((len(Temp),len(theta),len(phi)))

p1lndA1 = np.zeros((len(Temp),len(theta),len(phi)))
p2lndA1 = np.zeros((len(Temp),len(theta),len(phi)))

derP1A1 = np.zeros((len(Temp)-1,len(theta),len(phi)))
derP2A1 = np.zeros((len(Temp)-1,len(theta),len(phi)))

#################
p1dA2 = np.zeros((len(Temp),len(theta),len(phi)))
p2dA2 = np.zeros((len(Temp),len(theta),len(phi)))

p1lndA2 = np.zeros((len(Temp),len(theta),len(phi)))
p2lndA2 = np.zeros((len(Temp),len(theta),len(phi)))

derP1A2 = np.zeros((len(Temp)-1,len(theta),len(phi)))
derP2A2 = np.zeros((len(Temp)-1,len(theta),len(phi)))

#################
p1dA3 = np.zeros((len(Temp),len(theta),len(phi)))
p2dA3 = np.zeros((len(Temp),len(theta),len(phi)))

p1lndA3 = np.zeros((len(Temp),len(theta),len(phi)))
p2lndA3 = np.zeros((len(Temp),len(theta),len(phi)))

derP1A3 = np.zeros((len(Temp)-1,len(theta),len(phi)))
derP2A3 = np.zeros((len(Temp)-1,len(theta),len(phi)))

#################
p1dA4 = np.zeros((len(Temp),len(theta),len(phi)))
p2dA4 = np.zeros((len(Temp),len(theta),len(phi)))

p1lndA4 = np.zeros((len(Temp),len(theta),len(phi)))
p2lndA4 = np.zeros((len(Temp),len(theta),len(phi)))

derP1A4 = np.zeros((len(Temp)-1,len(theta),len(phi)))
derP2A4 = np.zeros((len(Temp)-1,len(theta),len(phi)))

#################
p1dA5 = np.zeros((len(Temp),len(theta),len(phi)))
p2dA5 = np.zeros((len(Temp),len(theta),len(phi)))

p1lndA5 = np.zeros((len(Temp),len(theta),len(phi)))
p2lndA5 = np.zeros((len(Temp),len(theta),len(phi)))

derP1A5 = np.zeros((len(Temp)-1,len(theta),len(phi)))
derP2A5 = np.zeros((len(Temp)-1,len(theta),len(phi)))

#################
F1A1 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A1 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA1 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA1 = np.zeros((len(Temp)-1,len(theta)))

QFIA1 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A2 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A2 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA2 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA2 = np.zeros((len(Temp)-1,len(theta)))

QFIA2 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A3 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A3 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA3 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA3 = np.zeros((len(Temp)-1,len(theta)))

QFIA3 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A4 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A4 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA4 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA4 = np.zeros((len(Temp)-1,len(theta)))

QFIA4 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A5 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A5 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA5 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA5 = np.zeros((len(Temp)-1,len(theta)))

QFIA5 = np.zeros((len(Temp)-1,len(range(0, n))))

QFI_Therm =  np.zeros((len(Temp)))

#################

PA = np.zeros(len(NN))
###############################################################################
#################################  Operadores  ################################

N=2

# System Alone:
sm= sigmap()
sp= sigmam()
sx= sigmax()
sy= sigmay()
sz= sigmaz()

#### System: 4 qubits
# qubit 1:
Sm1 = tensor(sigmap(),qeye(2),qeye(2),qeye(2),qeye(2))
Sp1 = tensor(sigmam(),qeye(2),qeye(2),qeye(2),qeye(2))
Sx1 = tensor(sigmax(),qeye(2),qeye(2),qeye(2),qeye(2))
Sy1 = tensor(sigmay(),qeye(2),qeye(2),qeye(2),qeye(2))
Sz1 = tensor(sigmaz(),qeye(2),qeye(2),qeye(2),qeye(2))

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

#################################### Cramer-Rao ###########################################
# TempMax = 0.3

# ntherm = 1/(exp((OmegaS)/(TempMax))-1)
# Fth = (1/(ntherm*(ntherm+1)*(2*ntherm+1)*(2*ntherm+1)))*((exp(OmegaS/TempMax)*OmegaS)/((-1+exp(OmegaS/TempMax))*(-1+exp(OmegaS/TempMax))*TempMax*TempMax))

###############################################################################
###################################  Hamiltonian ##############################

HP4A = -J * (Sz1* Sz2 + Sz2* Sz3+ Sz3* Sz4 + Sz4* Sz5)
HP4B = -g * (Sx1 + Sx2 + Sx3 + Sx4 + Sx5)

H = HP4A + HP4B


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
    
###############################################################################
###############################  Collapse Operator (Lindblad) ###########################

    C1S = np.sqrt(gamma*(nt+1))*Sm
    C2S = np.sqrt(gamma*nt)*Sp

    Clist = [C1S,C2S] 

###############################################################################
#################################  Initial State  #############################

    G = basis(2,0) # base: excited state
    E = basis(2,1) # base: ground state


#######  Ancilla 1
    A1_st = G
    
#######  Ancilla 2
    A2_st = G
    
#######  Ancilla 3
    A3_st = G
    
 #######  Ancilla 4   
    A4_st = G
    
 #######  Ancilla 5   
    A5_st = G
    
    psi = tensor(A1_st,A2_st,A3_st,A4_st,A5_st) # vetor

#######  Completo

    S0mat = psi*psi.dag() # matriz

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
###############################   POVM Ancilla 1   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1

            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()

            

            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()


            p1A1 = (P1a*Rhof1).tr() # probabilidade
            R1A1 = P1a*Rhof1



            p2A1 = (P2a*Rhof1).tr()
            R2A1 = P2a*Rhof1

        
            p1dA1[r,q,z] = real(p1A1)
            p2dA1[r,q,z] = real(p2A1)
        
            p1lnA1 = np.log(p1A1)
            p2lnA1 = np.log(p2A1)

            p1lndA1[r,q,z] = real(p1lnA1)
            p2lndA1[r,q,z] = real(p2lnA1)

###############################################################################                
###############################   POVM Ancilla 2   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1

            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()

        

            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()


            p1A2 = (P1a*Rhof2).tr()
            R1A2 = P1a*Rhof2

            p2A2 = (P2a*Rhof2).tr()
            R2A2 = P2a*Rhof2
            
        
            p1dA2[r,q,z] = real(p1A2)
            p2dA2[r,q,z] = real(p2A2)
        
            p1lnA2 = np.log(p1A2)
            p2lnA2 = np.log(p2A2)
            
            p1lndA2[r,q,z] = real(p1lnA2)
            p2lndA2[r,q,z] = real(p2lnA2)
                
###############################################################################                
###############################   POVM Ancilla 3   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1

            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()

        

            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()


            p1A3 = (P1a*Rhof3).tr()
            R1A3 = P1a*Rhof3

            p2A3 = (P2a*Rhof3).tr()
            R2A3 = P2a*Rhof3
           
        
            p1dA3[r,q,z] = real(p1A3)
            p2dA3[r,q,z] = real(p2A3)
        
            p1lnA3 = np.log(p1A3)
            p2lnA3 = np.log(p2A3)
            
            
            p1lndA3[r,q,z] = real(p1lnA3)
            p2lndA3[r,q,z] = real(p2lnA3)
                
###############################################################################                
###############################   POVM Ancilla 4   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1

            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()

        
            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()
        
            p1A4 = (P1a*Rhof4).tr()
            R1A4 = P1a*Rhof4
           
            p2A4 = (P2a*Rhof4).tr()
            R2A4 = P2a*Rhof4
            
            
            p1dA4[r,q,z] = real(p1A4)
            p2dA4[r,q,z] = real(p2A4)
        
            p1lnA4 = np.log(p1A4)
            p2lnA4 = np.log(p2A4)

            
            p1lndA4[r,q,z] = real(p1lnA4)
            p2lndA4[r,q,z] = real(p2lnA4)
            
###############################################################################                
###############################   POVM Ancilla 5   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1

            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()

        
            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()
        
            p1A5 = (P1a*Rhof5).tr()
            R1A5 = P1a*Rhof5
           
            p2A5 = (P2a*Rhof5).tr()
            R2A5 = P2a*Rhof5
            
            
            p1dA5[r,q,z] = real(p1A5)
            p2dA5[r,q,z] = real(p2A5)
        
            p1lnA5 = np.log(p1A5)
            p2lnA5 = np.log(p2A5)

            
            p1lndA5[r,q,z] = real(p1lnA5)
            p2lndA5[r,q,z] = real(p2lnA5)
                
###############################################################################
##################################   Ancilla 1  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):

        for p in range(len(phi)):
            derP1A1[r,i,p] = (p1lndA1[r+1,i,p] - p1lndA1[r,i,p])/dtSE
            derP2A1[r,i,p] = (p2lndA1[r+1,i,p] - p2lndA1[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 
    for i in range(len(theta)):

        for p in range(len(phi)):
            f1 = p1dA1[r,i,p]*(derP1A1[r,i,p])**2
            f2 = p2dA1[r,i,p]*(derP2A1[r,i,p])**2

            fx = f1+f2

    
            F1A1[r,i,p] = f1
            F2A1[r,i,p] = f2
    
            FxA1[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA1[r,i] = max(FxA1[r,i,:])

for r in range(len(Temp)-1):
        QFIA1[r] = max(QFI_tA1[r,:])
            

###############################################################################
##################################   Ancilla 2  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):

        for p in range(len(phi)):
            derP1A2[r,i,p] = (p1lndA2[r+1,i,p] - p1lndA2[r,i,p])/dtSE
            derP2A2[r,i,p] = (p2lndA2[r+1,i,p] - p2lndA2[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 
    for i in range(len(theta)):

        for p in range(len(phi)):
            f1 = p1dA2[r,i,p]*(derP1A2[r,i,p])**2
            f2 = p2dA2[r,i,p]*(derP2A2[r,i,p])**2
                    

            fx = f1+f2
            

    
            F1A2[r,i,p] = f1
            F2A2[r,i,p] = f2
    
            FxA2[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA2[r,i] = max(FxA2[r,i,:])
    
for r in range(len(Temp)-1):
    QFIA2[r] = max(QFI_tA2[r,:])
            
        
###############################################################################
##################################   Ancilla 3  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):

        for p in range(len(phi)):
            derP1A3[r,i,p] = (p1lndA3[r+1,i,p] - p1lndA3[r,i,p])/dtSE
            derP2A3[r,i,p] = (p2lndA3[r+1,i,p] - p2lndA3[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 
    for i in range(len(theta)):

        for p in range(len(phi)):
            f1 = p1dA3[r,i,p]*(derP1A3[r,i,p])**2
            f2 = p2dA3[r,i,p]*(derP2A3[r,i,p])**2
                    

            fx = f1+f2
            
    
            F1A3[r,i,p] = f1
            F2A3[r,i,p] = f2
    
            FxA3[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA3[r,i] = max(FxA3[r,i,:])

for r in range(len(Temp)-1):
    QFIA3[r] = max(QFI_tA3[r,:])


###############################################################################
##################################   Ancilla 4  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):

        for p in range(len(phi)):
            derP1A4[r,i,p] = (p1lndA4[r+1,i,p] - p1lndA4[r,i,p])/dtSE
            derP2A4[r,i,p] = (p2lndA4[r+1,i,p] - p2lndA4[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 

    for i in range(len(theta)):

        for p in range(len(phi)):
            f1 = p1dA4[r,i,p]*(derP1A4[r,i,p])**2
            f2 = p2dA4[r,i,p]*(derP2A4[r,i,p])**2
                    

            fx = f1+f2
            
    
            F1A4[r,i,p] = f1
            F2A4[r,i,p] = f2
    
            FxA4[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA4[r,i] = max(FxA4[r,i,:])

for r in range(len(Temp)-1):
    QFIA4[r] = max(QFI_tA4[r,:])


###############################################################################
##################################   Ancilla 5  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):

        for p in range(len(phi)):
            derP1A5[r,i,p] = (p1lndA5[r+1,i,p] - p1lndA5[r,i,p])/dtSE
            derP2A5[r,i,p] = (p2lndA5[r+1,i,p] - p2lndA5[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 

    for i in range(len(theta)):

        for p in range(len(phi)):
            f1 = p1dA5[r,i,p]*(derP1A5[r,i,p])**2
            f2 = p2dA5[r,i,p]*(derP2A5[r,i,p])**2
                    

            fx = f1+f2
            
    
            F1A5[r,i,p] = f1
            F2A5[r,i,p] = f2
    
            FxA5[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA5[r,i] = max(FxA5[r,i,:])

for r in range(len(Temp)-1):
    QFIA5[r] = max(QFI_tA5[r,:])



###############################################################################
FthSca = np.zeros((len(Temp)-1))
for r in range(len(Temp)-1):
    
    FthSca[r] = (1/(nthermo[r+1]*(nthermo[r+1]+1)*(2*nthermo[r+1]+1)**2))*((nthermo[r+1]-nthermo[r])/dTemp)**2

###############################################################################

Write_Outfile(ddTemp, QFIA1[:,-1], f'./Results/QFI_q1_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(ddTemp, QFIA2[:,-1], f'./Results/QFI_q2_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(ddTemp, QFIA3[:,-1], f'./Results/QFI_q3_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(ddTemp, QFIA4[:,-1], f'./Results/QFI_q4_g{g:.2f}_ttherm{tSEmax:.3f}.txt')
Write_Outfile(ddTemp, QFIA5[:,-1], f'./Results/QFI_q5_g{g:.2f}_ttherm{tSEmax:.3f}.txt')

###############################################################################
###############################################################################


tend = time.time()    # tempo final de processamento
delta = tend - tin    # funcao calculo do intervalo de tempo  de processamento
print (delta)


