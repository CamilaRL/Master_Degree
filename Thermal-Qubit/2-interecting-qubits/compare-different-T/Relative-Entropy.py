import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from qutip import *
import os


def Read_Density_Matrices(filename, dim):

    rhoTList = np.loadtxt(filename, dtype='complex')

    tlist = []
    rhot = []
    
    for t, rho in enumerate(rhoTList):
        
        tlist.append(t)
        rhot.append(Qobj(rho.reshape((dim, dim))))
    
    return tlist, rhot


def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))
    
    
def Entropia_Relativa_Populacoes(pi, pt):
    
    return pi*np.log(pi/pt) + (1-pi)*np.log((1-pi)/(1-pt))
    
    
def Derivada(xlist, ylist):
    
    ## https://www.youtube.com/watch?v=utRKIlOZbtw
    
    yprime = np.diff(ylist)/np.diff(xlist)
    xprime = []
    
    for i in range(len(yprime)):
        
        xtemp = (xlist[i+1] + xlist[i])/2
        xprime = np.append(xprime, xtemp)
    
    return xprime, yprime

    
  

######### MAIN #########

w0 = 1

Tf_qubit, Th, Tc = np.loadtxt('./DensityMatrices/temperature_qubit.txt', unpack=True, usecols=(1,2,3))

gList = [0, 0.8]
cList = ['min', 'max']

tlist = np.arange(0, 30, 0.01)

for i, g in enumerate(gList):
    
    p_final = pFunc(Tf_qubit[i], w0)
    
    for c in cList:
    
        Sc1 = []
        Sc2 = []
    
        tempo_index, rhof_q1 = Read_Density_Matrices(f'./DensityMatrices/rhof_q1_c{c}_g{g}.txt', 2)
        tempo_index, rhof_q2 = Read_Density_Matrices(f'./DensityMatrices/rhof_q2_c{c}_g{g}.txt', 2)
        
        for t in tempo_index:
        
            p1 = rhof_q1[t].full()[0][0]
            p2 = rhof_q2[t].full()[0][0]
            
            Sc1.append( Entropia_Relativa_Populacoes(p1, p_final) )
            Sc2.append( Entropia_Relativa_Populacoes(p2, p_final) )
        
        
        dtlist1, dSc1 = Derivada(tlist, Sc1)
        dtlist2, dSc2 = Derivada(tlist, Sc2)
        
        
        for j in range(len(dSc1)):
            
            dSc1[j] = - dSc1[j]
            dSc2[j] = - dSc2[j]
    
        ## save data    
    
        fSc1 = open(f'./Thermodynamics/Sr_q1_c{c}_g{g}.txt', 'w')
        fSc2 = open(f'./Thermodynamics/Sr_q2_c{c}_g{g}.txt', 'w')
        fdSc1 = open(f'./Thermodynamics/dSr_q1_c{c}_g{g}.txt', 'w')
        fdSc2 = open(f'./Thermodynamics/dSr_q2_c{c}_g{g}.txt', 'w')
        
        for t in range(len(tlist)):
            
            fSc1.write(f'{tlist[t]} {Sc1[t]}\n')
            fSc2.write(f'{tlist[t]} {Sc2[t]}\n')
            
        for t in range(len(dtlist1)):
            
            fdSc1.write(f'{dtlist1[t]} {dSc1[t]}\n')
            fdSc2.write(f'{dtlist2[t]} {dSc2[t]}\n')
    
        ## plot data
    
        plt.plot(tlist, Sc1, color='red', label='Qubit 1')
        plt.plot(tlist, Sc2, color='blue', label='Qubit 2')
        plt.ylabel('Relative Entropy')
        plt.xlabel('Time')
        plt.title(f'c{c} | g = {g}')
        plt.legend()
        plt.show()
        
        plt.plot(dtlist1, dSc1, color='red', label='Qubit 1')
        plt.plot(dtlist2, dSc2, color='blue', label='Qubit 2')
        plt.ylabel('Local Entropy Production Rate')
        plt.xlabel('Time')
        plt.title(f'c{c} | g = {g}')
        plt.legend()
        plt.show()

    
    
