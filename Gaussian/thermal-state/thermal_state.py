import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f


def Delta_Beta(w, beta_i, beta, gamma, t):

    fi = Distribution('bose', beta_i, w)
    
    ff = Distribution('bose', beta, w)
    
    delta_beta = (fi - ff) * np.exp(-gamma*t) + ff
    
    d_delta_beta = -gamma * (fi - ff) * np.exp(-gamma*t)
    
    return delta_beta, d_delta_beta
    

def Wigner_Fisher_Info(w, beta_i, beta_f, gamma, t):

    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    return 2 * (d_delta_beta/delta_beta)**2
    

def Velocity(Iw):

    return np.sqrt(Iw)
    

def Distance(w, beta_i, beta_f, gamma, t):

    fi = Distribution('bose', beta_i, w)
    
    ff = Distribution('bose', beta_f, w)

    L = np.sqrt(2) * abs(np.log(((fi - ff) * np.exp(-gamma*t) + ff) / fi))
    
    return L
    

def EntropyProduction(w, beta_i, beta_f, gamma, t):
    
    ff = Distribution('bose', beta_f, w)
    
    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    Sprod = - (d_delta_beta/ff) + (d_delta_beta/delta_beta)
    
    return Sprod


def ThermalKinematics(beta_i, beta_f, w, gamma, tlist):

    Iw = []
    Vw = []
    Lw = []
    Kevol = []
    Sprod = []

    for t in tlist:

        Iwt = Wigner_Fisher_Info(w, beta_i, beta_f, gamma, t)

        Vwt = Velocity(Iwt)
        
        Lwt = Distance(w, beta_i, beta_f, gamma, t)

        Iw.append(Iwt)
        Vw.append(Vwt)
        Lw.append(Lwt)

        Kevol.append(RelativeEntropy(beta_i, beta_f, w, gamma, t))
        
        Sprod.append(EntropyProduction(w, beta_i, beta_f, gamma, t))
        

    completion = []

    for i in range(len(Lw)):

        completion.append(Lw[i]/Lw[-1])
        
    return Iw, Vw, Lw, completion, Kevol, Sprod
    

def RelativeEntropy(beta_i, beta_f, w, gamma, t):

    delta_i, d_delta_i = Delta_Beta(w, beta_i, beta_f, gamma, t)
    ff = Distribution('bose', beta_f, w)
    
    return -1 + (delta_i/ff) + np.log(ff/delta_i)


def EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list):

    K_diff = lambda beta_i: RelativeEntropy(beta_i, beta_eq, w, gamma, 0) - Kinit

    K_diff_list = []
    for beta in beta_list:
    
        K_diff_list.append(K_diff(beta))

    ## hot
    for i in range(0, np.argmin(K_diff_list), 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_1 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K1 = RelativeEntropy(beta_1, beta_eq, w, gamma, 0)

    ## cold
    for i in range(np.argmin(K_diff_list), len(beta_list) - 1, 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_2 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K2 = RelativeEntropy(beta_2, beta_eq, w, gamma, 0)
    
    
    ## plot
    
    K = [RelativeEntropy(beta_i, beta_eq, w, gamma, 0) for beta_i in beta_list]
    
    plt.plot(beta_list, K, color='orange')
    plt.scatter([beta_eq], [RelativeEntropy(beta_eq, beta_eq, w, gamma, 0)], color='orange', label=f'Tw = {1/beta_eq:.3f}')
    plt.scatter([beta_2], [K2], color='blue', label=f'Tc = {1/beta_2:.3f}')
    plt.scatter([beta_1], [K1], color='red', label=f'Th = {1/beta_1:.3f}')
    plt.hlines(y=Kinit, xmin=min(beta_list), xmax=max(beta_list), color='black', label='Initial Relative Entropy')
    plt.xlabel(r'$\beta$')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()
    
    
    return beta_1, 1/beta_1, K1, beta_2, 1/beta_2, K2
        
        





########### MAIN ###########


## Parameters

w = 1
gamma = 0.1

Kinit = 1

Teq = 2
beta_eq = 1/Teq

beta_list = np.arange(0.1, 5, 0.01)

beta_hot, Thot, Khot, beta_cold, Tcold, Kcold = EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list)

tlist = np.arange(0, 80, 0.1)


## Thermal Kinematics

Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating = ThermalKinematics(beta_cold, beta_eq, w, gamma, tlist)
Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling = ThermalKinematics(beta_hot, beta_eq, w, gamma, tlist)
    
    
## Plots

### Wigner Fisher Information

plt.figure(figsize=(10,5))

plt.subplot(121)
plt.plot(tlist, Iw_heating, color='red', label='Heating')
plt.ylabel('Wigner Fisher Information', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.subplot(122)
plt.plot(tlist, Iw_cooling, color='blue', label='Cooling')
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.tight_layout()
plt.show()

### Velocity

plt.figure(figsize=(10,5))

plt.subplot(121)
plt.plot(tlist, Vw_heating, color='red', label='Heating')
plt.ylabel('Statistical Velocity', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.subplot(122)
plt.plot(tlist, Vw_cooling, color='blue', label='Cooling')
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.tight_layout()
plt.show()

### Distance

plt.plot(tlist, Lw_heating, color='red', label='Heating')
plt.plot(tlist, Lw_cooling, color='blue', label='Cooling')
plt.ylabel('Statistical Distance', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.legend(fontsize=12)
plt.show()

### Completion

plt.plot(tlist, completion_heating, color='red', label='Heating')
plt.plot(tlist, completion_cooling, color='blue', label='Cooling')
plt.ylabel('Degree of Completion', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.legend(fontsize=12)
plt.show()

### Relative Entropy

plt.plot(tlist, Kevol_heating, color='red', label='Heating')
plt.plot(tlist, Kevol_cooling, color='blue', label='Cooling')
plt.ylabel('Relative Entropy', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.legend(fontsize=12)
plt.show()

### Entropy Production

plt.plot(tlist, Sprod_heating, color='red', label='Heating')
plt.plot(tlist, Sprod_cooling, color='blue', label='Cooling')
plt.ylabel('Entropy Production', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.legend(fontsize=12)
plt.show()








