import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import math
import os



def pFunc(T, w0):

    return np.exp(w0/(2*T))/(2*np.cosh(w0/(2*T)))
    
    
def Relacao_Igualdade_Srelativa(Sr_init, w0, Tf, dT):

    lado_esquerdo = Sr_init - np.log(1 + np.exp(w0/Tf))
    
    func_direito = lambda T0: ((w0/T0) - (w0/Tf))/(1 + np.exp(-w0/T0)) - np.log(1 + np.exp(w0/T0)) - lado_esquerdo
    
    Tlist = np.arange(-Tf-50, Tf+50, dT)
    
    lado_direito = func_direito(Tlist)
    
    y_cruzamento = []
    T_cruzamento = []
    
    for i in range(0, len(Tlist)-1, 1):
    
        if lado_direito[i] * lado_direito[i+1] < 0:
            
            Ti = brentq(func_direito, Tlist[i], Tlist[i+1])
            
            y_cruzamento.append(func_direito(Ti))
            T_cruzamento.append(Ti)
            
    print(Tf)
    print(T_cruzamento)
    
    plt.plot(Tlist, lado_direito)
    plt.scatter(T_cruzamento, y_cruzamento)
    plt.show()


def Interseccao_Temperatura_Inicial(w0, p, Tlist):
    
    find_T = lambda T: pFunc(T, w0) - p
        
    Ti = brentq(find_T, Tlist[0], Tlist[-1])  # busca em todo o range de T
    
    return Ti


##### MAIN #####

## parametros

w0 = 2

gamma = 1

p_final = 0.4

Sr_init = 0.1 ### se alterar, deve mudar os ranges de temperatura 

tlist = np.arange(0, 10, 0.01)


## estados equidistantes

print('Temperatura e Populações')

T_list = np.arange(-50, -0.0001, 0.0001)

Tw = Interseccao_Temperatura_Inicial(w0, p_final, T_list)

Relacao_Igualdade_Srelativa(Sr_init, w0, Tw, dT=0.0001)







