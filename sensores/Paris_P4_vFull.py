from qutip import *
from math import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import time
from numpy import linalg as LA
from numpy.linalg import norm
import scipy.constants as constant
# import matplotlib.pyplot as plt

tin = time.time()

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

g = 1.0  # Sistema-Campo
J = 1.0       # Sistema-sistema

# OmegaS = g#1E10*(T/2)  # System: energy level splitting

# raz_To= (constant.k*T)/(constant.h*Omega)
# raz_To= (kb*T)/((hp/2*pi)*Omega)

# print(raz_To)

tT = 100
TempMin = 0.001
TempMax = 2.0
Temp = linspace(TempMin,TempMax,tT) # Temperature

dtT = tT - 1
ddTemp = linspace(TempMin,TempMax,dtT) # Temperature

dTemp = Temp[1] - Temp[0]

tp = 100 # Step
tSA1 = linspace(0.00001,pi/2,tp)  # Time S-A
# tSA1 = linspace(0.00001,pi/100,tp)  # Time S-A1
tA1A2 = linspace(0.00001,pi/2,tp)  # Time A1-A2
tA2A3 = linspace(0.00001,pi/2,tp)  # Time A2-A3
tA3A4 = linspace(0.00001,pi/2,tp)  # Time A3-A4
tA4A5 = linspace(0.00001,pi/2,tp)  # Time A4-A5
tA5A6 = linspace(0.00001,pi/2,tp)  # Time A5-A6


# tSE = linspace(0.00001,0.1,tp)  # Time S-E
tSE = linspace(0.00001,pi/2,tp)  # Time S-E

td = linspace(0.0,5.0,tp-1)
# dt = 5/tp
dtSA1 = tSA1[1]-tSA1[0]
dtA1A2 = tA1A2[1]-tA1A2[0]

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

Coh_S1 = np.zeros((len(Temp),len(range(0, n))))
Coh_S2 = np.zeros((len(Temp),len(range(0, n))))
Coh_S3 = np.zeros((len(Temp),len(range(0, n))))
Coh_S4 = np.zeros((len(Temp),len(range(0, n))))

Coh_A1U1 = np.zeros((len(Temp),len(range(0, n))))
Coh_A2U2 = np.zeros((len(Temp),len(range(0, n))))
Coh_A3U2 = np.zeros((len(Temp),len(range(0, n))))
Coh_A4U2 = np.zeros((len(Temp),len(range(0, n))))

ES = np.zeros((len(Temp),len(range(0, n))))


nthermo = np.zeros((len(Temp)))

FthGab = np.zeros((len(Temp)))

EntR1_A1 = np.zeros((len(Temp),len(theta),len(phi)))
EntR2_A1 = np.zeros((len(Temp),len(theta),len(phi)))

EntR1_A2 = np.zeros((len(Temp),len(theta),len(phi)))
EntR2_A2 = np.zeros((len(Temp),len(theta),len(phi)))

EntR1_A3 = np.zeros((len(Temp),len(theta),len(phi)))
EntR2_A3 = np.zeros((len(Temp),len(theta),len(phi)))

EntR1_A4 = np.zeros((len(Temp),len(theta),len(phi)))
EntR2_A4 = np.zeros((len(Temp),len(theta),len(phi)))


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
F1A1 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A1 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA1 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA1 = np.zeros((len(Temp)-1,len(theta)))
# QFI_t = np.zeros((len(Temp)-1,len(phi)))

QFIA1 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A2 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A2 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA2 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA2 = np.zeros((len(Temp)-1,len(theta)))
# QFI_t = np.zeros((len(Temp)-1,len(phi)))

QFIA2 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A3 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A3 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA3 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA3 = np.zeros((len(Temp)-1,len(theta)))
# QFI_t = np.zeros((len(Temp)-1,len(phi)))

QFIA3 = np.zeros((len(Temp)-1,len(range(0, n))))

#################
F1A4 =  np.zeros((len(Temp)-1,len(theta),len(phi)))
F2A4 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

FxA4 =  np.zeros((len(Temp)-1,len(theta),len(phi)))

QFI_tA4 = np.zeros((len(Temp)-1,len(theta)))
# QFI_t = np.zeros((len(Temp)-1,len(phi)))

QFIA4 = np.zeros((len(Temp)-1,len(range(0, n))))

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
# TempMax = 0.3

# ntherm = 1/(exp((OmegaS)/(TempMax))-1)
# Fth = (1/(ntherm*(ntherm+1)*(2*ntherm+1)*(2*ntherm+1)))*((exp(OmegaS/TempMax)*OmegaS)/((-1+exp(OmegaS/TempMax))*(-1+exp(OmegaS/TempMax))*TempMax*TempMax))

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
r = -1
# q = 0
for T in Temp:
    # print(T)
    
    # q = q + 1
    r = r + 1
    print(r)
###############################################################################
###############################  Thermal Number  ##############################
    # nt = 1/(exp((hp*Omega)/(kb*T))-1)
    nt = 1/(exp((OmegaS)/(T))-1)
    nthermo[r] = nt
    # print(nt)
    # sech = 1/cosh(x)
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
    # A1_st = 1/(sqrt(2))*(G+E)
    A1_mat = A1_st*A1_st.dag()
    # print(A_mat)
#######  Ancilla 2
    A2_st = G
    # A2_st = 1/(sqrt(2))*(G+E)
    A2_mat = A2_st*A2_st.dag()
    # # print(A_mat)
#######  Ancilla 3
    A3_st = G
    # A2_st = 1/(sqrt(2))*(G+E)
    A3_mat = A3_st*A3_st.dag()
    # # print(A_mat)
 #######  Ancilla 4   
    A4_st = G
    # A2_st = 1/(sqrt(2))*(G+E)
    A4_mat = A4_st*A4_st.dag()
    # # print(A_mat)

    
    psi = tensor(A1_st,A2_st,A3_st,A4_st)
    # rho0 = psi*psi.dag()
#######  Completo
    # psi = tensor(S_st[-1],A1_st,A2_st)
    S0mat = psi*psi.dag()

    medataSE = mesolve(H,S0mat,tSE,[Clist],[]) # Master equation evolution - Sistema + Ambiente
    rho0 = medataSE.states # Take matrices in each time
    medataH1 = mesolve(H,S0mat,tSE,[Clist],[Sz1]) # Master equation evolution - Sistema + Ambiente
    Exp_rho0 = medataH1.expect[0] # Take matrices in each time
    # print (Exp_rho0)
    
    rhof = rho0[-1] # Take matrices in each time
    Rhof1 = rho0[-1].ptrace([0])
    Rhof2 = rho0[-1].ptrace([1])
    Rhof3 = rho0[-1].ptrace([2])
    Rhof4 = rho0[-1].ptrace([3])
    # print(A_evo)
    
    medataH2 = mesolve(H,S0mat,tSE,[Clist],[Sz1*Sz1]) # Master equation evolution - Sistema + Ambiente
    Exp_rho02 = medataH2.expect[0]
    # print (Exp_rho02)
    
    # QFI_Therm[r] = (Exp_rho02[-1] - (Exp_rho0[-1]*Exp_rho0[-1]))/(T**4)
    Therm1 = (Exp_rho02[-1] - (Exp_rho0[-1]*Exp_rho0[-1]))
    # print (Therm1)
    # QFI_Therm[r] = np.divide(Therm1,T**4)
    QFI_Therm[r] = Therm1/(T**4)
    
    


    # rho0 = tensor(exptSE[-1],A1_mat,A2_mat,A3_mat,A4_mat,A5_mat,A6_mat)
        
###############################################################################    
###############################   POVM   ############################  
    q = -1
    for i in range(len(theta)):
        q = q + 1
        z = -1
        for p in range(len(phi)):
            z = z + 1
            # mat1 = Qobj([[exp(-1j*phi[p])], [exp(1j*phi[p])*cos(theta[i])]])
            mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
            P1a = mat1*mat1.dag()
            P1 = tensor(P1a,P1a,P1a,P1a)
            
            # P2 = -1j*exp(1j*phi[p])*sin(theta[i])*sm
            mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
            P2a = mat2*mat2.dag()
            P2 = tensor(P2a,P2a,P2a,P2a)

            # p1A1 = (P1a*Rhof1).tr()
            p1A = (P1*rhof).tr()
            # R1A1 = P1a*Rhof1
            R1A = P1*rhof
            # p1A1 = (P1a*Rho_A1U2).tr()
            # p1A1 = (P1a*Rho_A1ciclo).tr()
            # print(p1)
            p2A = (P2*rhof).tr()
            R2A = P2*rhof
            # p2A1 = (P2a*Rho_A1U2).tr()
            # p2A1 = (P2a*Rho_A1ciclo).tr()
            
            # EntR1_A1[r,s,q,z] = entropy_relative(Rho_A1U1, R1A1, base=e, sparse=False)
            # EntR2_A1[r,s,q,z] = entropy_relative(Rho_A1U1, R2A1, base=e, sparse=False)
        
#             # p1 = (P1*exptSA[-1]).tr()
#             # # print(p1)
#             # p2 = (P2*exptSA[-1]).tr()
        
            p1dA1[r,q,z] = real(p1A)
            p2dA1[r,q,z] = real(p2A)
        
            p1lnA1 = np.log(p1A)
            p2lnA1 = np.log(p2A)
            # print(p1ln)
            
            p1lndA1[r,q,z] = p1lnA1
            p2lndA1[r,q,z] = p2lnA1

# ###############################################################################                
# ###############################   POVM Ancilla 2   ############################  
#     q = -1
#     for i in range(len(theta)):
#         q = q + 1
#         z = -1
#         for p in range(len(phi)):
#             z = z + 1
#             # mat1 = Qobj([[exp(-1j*phi[p])], [exp(1j*phi[p])*cos(theta[i])]])
#             mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
#             P1a = mat1*mat1.dag()
#         # P1 = tensor(qeye(2),P1a)
        
#         # P2 = -1j*exp(1j*phi[p])*sin(theta[i])*sm
#             mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
#             P2a = mat2*mat2.dag()
#         # P2 = tensor(qeye(2),P2a)

#             p1A2 = (P1a*Rhof2).tr()
#             R1A2 = P1a*Rhof2
#             # print(p1)
#             p2A2 = (P2a*Rhof2).tr()
#             R2A2 = P2a*Rhof2
            
#             # EntR1_A2[r,q,z] = entropy_relative(Rho_A2U2, R1A2, base=e, sparse=False)
#             # EntR2_A2[r,q,z] = entropy_relative(Rho_A2U2, R2A2, base=e, sparse=False)
        
# #             # p1 = (P1*exptSA[-1]).tr()
# #             # # print(p1)
# #             # p2 = (P2*exptSA[-1]).tr()
        
#             p1dA2[r,q,z] = real(p1A2)
#             p2dA2[r,q,z] = real(p2A2)
        
#             p1lnA2 = np.log(p1A2)
#             p2lnA2 = np.log(p2A2)
#             # print(p1ln)
            
#             p1lndA2[r,q,z] = p1lnA2
#             p2lndA2[r,q,z] = p2lnA2
                
# ###############################################################################                
# ###############################   POVM Ancilla 3   ############################  
#     q = -1
#     for i in range(len(theta)):
#         q = q + 1
#         z = -1
#         for p in range(len(phi)):
#             z = z + 1
#             # mat1 = Qobj([[exp(-1j*phi[p])], [exp(1j*phi[p])*cos(theta[i])]])
#             mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
#             P1a = mat1*mat1.dag()
#         # P1 = tensor(qeye(2),P1a)
        
#         # P2 = -1j*exp(1j*phi[p])*sin(theta[i])*sm
#             mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
#             P2a = mat2*mat2.dag()
#         # P2 = tensor(qeye(2),P2a)

#             p1A3 = (P1a*Rhof3).tr()
#             R1A3 = P1a*Rhof3
#             # print(p1)
#             p2A3 = (P2a*Rhof3).tr()
#             R2A3 = P2a*Rhof3
            
#             # EntR1_A3[r,s,q,z] = entropy_relative(Rho_A3U2, R1A3, base=e, sparse=False)
#             # EntR2_A3[r,s,q,z] = entropy_relative(Rho_A3U2, R2A3, base=e, sparse=False)
        
# #             # p1 = (P1*exptSA[-1]).tr()
# #             # # print(p1)
# #             # p2 = (P2*exptSA[-1]).tr()
        
#             p1dA3[r,q,z] = real(p1A3)
#             p2dA3[r,q,z] = real(p2A3)
        
#             p1lnA3 = np.log(p1A3)
#             p2lnA3 = np.log(p2A3)
#             # print(p1ln)
            
#             p1lndA3[r,q,z] = p1lnA3
#             p2lndA3[r,q,z] = p2lnA3
                
# ###############################################################################                
# ###############################   POVM Ancilla 4   ############################  
#     q = -1
#     for i in range(len(theta)):
#         q = q + 1
#         z = -1
#         for p in range(len(phi)):
#             z = z + 1
#             # mat1 = Qobj([[exp(-1j*phi[p])], [exp(1j*phi[p])*cos(theta[i])]])
#             mat1 = Qobj([[cos(theta[i])], [exp(1j*phi[p])*sin(theta[i])]])
#             P1a = mat1*mat1.dag()
            
#         # P1 = tensor(qeye(2),P1a)
        
#         # P2 = -1j*exp(1j*phi[p])*sin(theta[i])*sm
#             mat2 = Qobj([[exp(-1j*phi[p])*sin(theta[i])],[-cos(theta[i])]])
#             P2a = mat2*mat2.dag()
        
#             p1A4 = (P1a*Rhof4).tr()
#             R1A4 = P1a*Rhof4
#             # print(p1)
#             p2A4 = (P2a*Rhof4).tr()
#             R2A4 = P2a*Rhof4
            
#             # EntR1_A4[r,s,q,z] = entropy_relative(Rho_A4U2, R1A4, base=e, sparse=False)
#             # EntR2_A4[r,s,q,z] = entropy_relative(Rho_A4U2, R2A4, base=e, sparse=False)
                    
#             p1dA4[r,q,z] = real(p1A4)
#             p2dA4[r,q,z] = real(p2A4)
        
#             p1lnA4 = np.log(p1A4)
            # p2lnA4 = np.log(p2A4)
            # # print(p1ln)
            
            # p1lndA4[r,q,z] = p1lnA4
            # p2lndA4[r,q,z] = p2lnA4
                
###############################################################################
##################################   Ancilla 1  ###############################
for r in range(len(Temp)-1):
    for i in range(len(theta)):
                # print(i)
        for p in range(len(phi)):
            derP1A1[r,i,p] = (p1lndA1[r+1,i,p] - p1lndA1[r,i,p])/dtSE
            derP2A1[r,i,p] = (p2lndA1[r+1,i,p] - p2lndA1[r,i,p])/dtSE
    
    
for r in range(len(Temp)-1): 
    for i in range(len(theta)):
                # print(i)
        for p in range(len(phi)):
            f1 = p1dA1[r,i,p]*(derP1A1[r,i,p])**2
            f2 = p2dA1[r,i,p]*(derP2A1[r,i,p])**2
                    # print(f1)
            fx = f1+f2
                    # FI = max(fx)
    
            F1A1[r,i,p] = f1
            F2A1[r,i,p] = f2
    
            FxA1[r,i,p] = fx
    
    
for r in range(len(Temp)-1):
    for i in range(len(theta)):
        QFI_tA1[r,i] = max(FxA1[r,i,:])

for r in range(len(Temp)-1):
        QFIA1[r] = max(QFI_tA1[r,:])
            

# plot(NN,QFIA1[-1,:])
# plot(ddTemp,QFIA1[:,-1])

# ###############################################################################
# ##################################   Ancilla 2  ###############################
# for r in range(len(Temp)-1):
#     for i in range(len(theta)):
#                 # print(i)
#         for p in range(len(phi)):
#             derP1A2[r,i,p] = (p1lndA2[r+1,i,p] - p1lndA2[r,i,p])/dtSE
#             derP2A2[r,i,p] = (p2lndA2[r+1,i,p] - p2lndA2[r,i,p])/dtSE
    
    
# for r in range(len(Temp)-1): 
# # print(r)
#     for i in range(len(theta)):
#                 # print(i)
#         for p in range(len(phi)):
#             f1 = p1dA2[r,i,p]*(derP1A2[r,i,p])**2
#             f2 = p2dA2[r,i,p]*(derP2A2[r,i,p])**2
                    
#                     # print(f1)
#             fx = f1+f2
            
#                     # FI = max(fx)
    
#             F1A2[r,i,p] = f1
#             F2A2[r,i,p] = f2
    
#             FxA2[r,i,p] = fx
    
    
# for r in range(len(Temp)-1):
#     for i in range(len(theta)):
#         QFI_tA2[r,i] = max(FxA2[r,i,:])
    
# for r in range(len(Temp)-1):
#     QFIA2[r] = max(QFI_tA2[r,:])
            
        
# ###############################################################################
# ##################################   Ancilla 3  ###############################
# for r in range(len(Temp)-1):
#     for i in range(len(theta)):
#                 # print(i)
#         for p in range(len(phi)):
#             derP1A3[r,i,p] = (p1lndA3[r+1,i,p] - p1lndA3[r,i,p])/dtSE
#             derP2A3[r,i,p] = (p2lndA3[r+1,i,p] - p2lndA3[r,i,p])/dtSE
    
    
# for r in range(len(Temp)-1): 
#     for i in range(len(theta)):
#                 # print(i)
#         for p in range(len(phi)):
#             f1 = p1dA3[r,i,p]*(derP1A3[r,i,p])**2
#             f2 = p2dA3[r,i,p]*(derP2A3[r,i,p])**2
                    
#                     # print(f1)
#             fx = f1+f2
            
#                     # FI = max(fx)
    
#             F1A3[r,i,p] = f1
#             F2A3[r,i,p] = f2
    
#             FxA3[r,i,p] = fx
    
    
# for r in range(len(Temp)-1):
#     for i in range(len(theta)):
#         QFI_tA3[r,i] = max(FxA3[r,i,:])

# for r in range(len(Temp)-1):
#     QFIA3[r] = max(QFI_tA3[r,:])

# # plot(NN,QFIA2[-1,:],'--r')
# # plot(ddTemp,QFIA2[:,-1],'--r')


# ###############################################################################
# ##################################   Ancilla 4  ###############################
# for r in range(len(Temp)-1):
#     for i in range(len(theta)):
#                 # print(i)
#         for p in range(len(phi)):
#             derP1A4[r,i,p] = (p1lndA4[r+1,i,p] - p1lndA4[r,i,p])/dtSE
#             derP2A4[r,i,p] = (p2lndA4[r+1,i,p] - p2lndA4[r,i,p])/dtSE
    
    
# for r in range(len(Temp)-1): 
# # print(r)
#     for i in range(len(theta)):
#                 # print(i)
#         for p in range(len(phi)):
#             f1 = p1dA4[r,i,p]*(derP1A4[r,i,p])**2
#             f2 = p2dA4[r,i,p]*(derP2A4[r,i,p])**2
                    
#                     # print(f1)
#             fx = f1+f2
            
#                     # FI = max(fx)
    
#             F1A4[r,i,p] = f1
#             F2A4[r,i,p] = f2
    
#             FxA4[r,i,p] = fx
    
    
# for r in range(len(Temp)-1):
#     for i in range(len(theta)):
#         QFI_tA4[r,i] = max(FxA4[r,i,:])

# for r in range(len(Temp)-1):
#     QFIA4[r] = max(QFI_tA4[r,:])

# plot(NN,QFIA2[-1,:],'--r')
# plot(ddTemp,QFIA2[:,-1],'--r')


###############################################################################
FthSca = np.zeros((len(Temp)-1))
for r in range(len(Temp)-1):
    
    FthSca[r] = (1/(nthermo[r+1]*(nthermo[r+1]+1)*(2*nthermo[r+1]+1)**2))*((nthermo[r+1]-nthermo[r])/dTemp)**2

###############################################################################
plt.plot(ddTemp,QFIA1[:,-1])#,'b',ddTemp,QFIA2[:,-1],'--r',ddTemp,QFIA3[:,-1],'-.k',ddTemp,QFIA4[:,-1],':g')
# plot(NN,QFIA1[12,:],NN,QFIA2[12,:],'--r',NN,QFIA3[12,:],'-.k',NN,QFIA4[12,:],':m',NN,QFIA5[12,:],NN,QFIA6[12,:],'--g')
plt.plot()
# nCam = 4
# ACam = linspace(1,nCam,nCam)  # Time S-A

# # MaxQFI1 = max(QFIA1[:,-1])
# # MaxQFI2 = max(QFIA2[:,-1])
# # MaxQFI3 = max(QFIA3[:,-1])
# # MaxQFI4 = max(QFIA4[:,-1])
# # MaxQFI5 = max(QFIA5[:,-1])

# MaxQFI1 = max(QFIA1[12,:])
# MaxQFI2 = max(QFIA2[12,:])
# MaxQFI3 = max(QFIA3[12,:])
# MaxQFI4 = max(QFIA4[12,:])

# QFImax = np.zeros(len(ACam))

# QFImax[0] = MaxQFI1
# QFImax[1] = MaxQFI2
# QFImax[2] = MaxQFI3
# QFImax[3] = MaxQFI4

# plot(ACam,QFImax)

# plot(ddTemp,QFIA1[:,-1]/FthSca,ddTemp,QFIA2[:,-1]/FthSca,'--r')

# plot(NN,QFIA1[13,:]/FthSca[13],NN,QFIA2[13,:]/FthSca[13],'--r')
# plot(NN,QFIA1[13,:]/(FthGab[13]/n),NN,QFIA2[13,:]/(FthGab[13]/n),'--r')
# plot(NN,QFIA1[13,:]/FthGab[13],NN,QFIA2[13,:]/FthGab[13],'--r')
# plot(NN,QFIA1[13,:]/(FthSca[13]/n),NN,QFIA2[13,:]/(FthSca[13]/n),'--r')

# plt.xscale("log")
# plt.yscale("log")
# plt.plot(NN,QFIA1[13,:]/(FthGab[13]/n),NN,QFIA2[13,:]/(FthGab[13]/n),'--r')
# plt.plot(NN,QFIA1[12,:]/FthGab[12],NN,QFIA2[12,:]/FthGab[12],'--r')
# plt.plot(NN,QFIA1[13,:]/FthGab[13],NN,QFIA2[13,:]/FthGab[13],'--r')

# plt.plot(NN,QFIA1[13,:]/FthSca[13],NN,QFIA2[13,:]/FthSca[13],'--r')

###############################################################################
###############################################################################


tend = time.time()    # tempo final de processamento
delta = tend - tin    # funcao calculo do intervalo de tempo  de processamento
print (delta)


