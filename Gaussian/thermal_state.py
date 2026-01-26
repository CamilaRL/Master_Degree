import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def f(w0, beta):
    
    return 1/(np.exp(w0 * beta) - 1)


def Statistical_Velocity(w0, beta_i, beta_f, gamma, t):

    return np.sqrt(2) * (-gamma * np.exp(-gamma*t) * (f(w0, beta_i) - f(w0, beta_f)))/(np.exp(-gamma*t) * (f(w0, beta_i) - f(w0, beta_f)) - f(w0, beta_f))
    
    
def Statistical_Length(w0, beta_i, beta_f, gamma, t):
    
    return np.sqrt(2) * abs(np.log((np.exp(-gamma*t) * (f(w0, beta_i) - f(w0, beta_f)) + f(w0, beta_f))/f(w0, beta_i)))


def Degree_of_Completion(Lt, Lf):
    
    return Lt/Lf


def ThermalKinematics(w0, gamma, tlist, beta_i, beta_f):
    
    vlist = []
    Llist = []

    for t in tlist:
        
        vt = Statistical_Velocity(w0, beta_i, beta_f, gamma, t)
        
        Lt = Statistical_Length(w0, beta_i, beta_f, gamma, t)

        vlist.append(vt)
        Llist.append(Lt)

    philist = []

    for Lt in Llist:
        
        philist.append(Degree_of_Completion(Lt, Llist[-1]))

    plt.plot(tlist, vlist)
    plt.ylabel('Velocity')
    plt.xlabel('Time')
    plt.show()
    
    plt.plot(tlist, Llist)
    plt.ylabel('Position')
    plt.xlabel('Time')
    plt.show()
    
    plt.plot(tlist, philist)
    plt.ylabel('Degree of Completion')
    plt.xlabel('Time')
    plt.show()
    

def Entropia_Relativa(w0, beta_i, beta_f):

    delta_beta_i = f(w0, beta_i)
    delta_beta_f = f(w0, beta_f)

    return -1 + 0.5*(2*delta_beta_i/delta_beta_f + np.log((delta_beta_f**2)/(delta_beta_i**2)))


def Equidistant_States(w0, beta_f, beta_list, Sr_init):
    
    f_diff =  lambda beta: Entropia_Relativa(w0, beta, beta_f) - Sr_init
    
    Sr_diff = f_diff(np.array(beta_list))
    
    ## hot
    for i in range(0, np.argmin(Sr_diff), 1):

        if Sr_diff[i] * Sr_diff[i+1] < 0:  # houve cruzamento

            beta_1 = brentq(f_diff, beta_list[i], beta_list[i+1])  # raiz exata
            Sr1 = Entropia_Relativa(w0, beta_1, beta_f)

    ## cold
    for i in range(np.argmin(Sr_diff), len(beta_list) - 1, 1):

        if Sr_diff[i] * Sr_diff[i+1] < 0:  # houve cruzamento

            beta_2 = brentq(f_diff, beta_list[i], beta_list[i+1])  # raiz exata
            Sr2 = Entropia_Relativa(w0, beta_2, beta_f)
    
    
    return beta_1, 1/beta_1, Sr1, beta_2, 1/beta_2, Sr2






### MAIN ###

w0 = 1
gamma = 0.1
Tf = 2
beta_f = 1/Tf
Sr_init = 0.5

tlist = np.arange(0, 20, 0.1)


## Equidistant

beta_list = np.arange(0.2, 2, 0.01)

Sr = [Entropia_Relativa(w0, beta_i, beta_f) for beta_i in beta_list]

beta_h, Th, Srh, beta_c, Tc, Src = Equidistant_States(w0, beta_f, beta_list, Sr_init)


plt.plot(beta_list, Sr, color='orange')
plt.scatter([beta_f], [0], color='orange', label=f'Tw = {Tf:.3f}')
plt.scatter([beta_c], [Src], color='blue', label=f'Tc = {Tc:.3f}')
plt.scatter([beta_h], [Srh], color='red', label=f'Th = {Th:.3f}')
plt.hlines(y=Sr_init, xmin=min(beta_list), xmax=max(beta_list), color='black', label='Initial Relative Entropy')
plt.xlabel(r'$\beta$')
plt.ylabel('Relative Entropy')
plt.legend()
plt.show()


## Heating
ThermalKinematics(w0, gamma, tlist, beta_c, beta_f)

## Cooling
ThermalKinematics(w0, gamma, tlist, beta_h, beta_f)
