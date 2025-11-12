import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import os


def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))


def Interseccao_Temperatura(p):
    
    find_T = lambda T: pFunc(T, w0=2) - p
        
    T = brentq(find_T, 0.01, 50)  # busca em todo o range de T
    
    return T


####### MAIN #######


modo = 'Cooling'

cmod = np.loadtxt(f'./DensityMatrices/cmod_{modo}.txt', unpack=True)

cmod_extremes = [min(cmod), max(cmod)]

tlist = np.arange(0.005, 2, 0.001)


for c in cmod_extremes:

    Tt = []
    coherence = []

    path = f'./DensityMatrices/{modo}/c_{c}/'

    rhoname = 'rhof_q1_t'

    for arquivo in os.listdir(path):
    
        if rhoname in arquivo:
            
            rhot = np.loadtxt(path + arquivo, unpack=True, dtype='complex')
        
            Tt.append( Interseccao_Temperatura(rhot[0][0].real) )
            
            coherence.append( abs(rhot[0][1])/2 )
            
            
    plt.scatter(tlist, Tt, s=2)
    plt.title(f'|c| = {c}')
    plt.ylabel('Temperature')
    plt.xlabel('Time')
    plt.show()

    plt.plot(tlist, coherence)
    plt.title(f'|c| = {c}')
    plt.ylabel('Coherence')
    plt.xlabel('Time')
    plt.show()
