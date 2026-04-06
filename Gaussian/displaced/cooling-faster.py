import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import math


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
    

def Wigner_Fisher_Info(mu, w, beta_i, beta_f, gamma, t):

    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    return 2 * (d_delta_beta/delta_beta)**2 + 2 * (w**2 + (gamma**2)/4) * (abs(mu)**2)*np.exp(-gamma*t)/delta_beta
    

def Velocity(Iw):

    return np.sqrt(Iw)
    

def Distance(vList, dt):
    
    L = np.cumsum(vList) * dt
    
    return L
    

def EntropyProduction(mu, w, beta_i, beta_f, gamma, t):
    
    ff = Distribution('bose', beta_f, w)
    
    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    Sprod = (gamma*(abs(mu)**2)*np.exp(-gamma*t) - d_delta_beta)/ff + (d_delta_beta/delta_beta)
    
    return Sprod


def ThermalKinematics(mu, beta_i, beta_f, w, gamma, tlist):

    Iw = []
    Vw = []
    Kevol = []
    Sprod = []

    for t in tlist:

        Iwt = Wigner_Fisher_Info(mu, w, beta_i, beta_f, gamma, t)

        Vwt = Velocity(Iwt)

        Iw.append(Iwt)
        Vw.append(Vwt)

        Kevol.append(RelativeEntropy(mu, beta_i, beta_f, w, gamma, t))
        
        Sprod.append(EntropyProduction(mu, w, beta_i, beta_f, gamma, t))
    
    Lw = Distance(Vw, tlist[1]-tlist[0]) 

    completion = []

    for i in range(len(Lw)):

        completion.append(Lw[i]/Lw[-1])
        
    return Iw, Vw, Lw, completion, Kevol, Sprod
    

def RelativeEntropy(mu, beta_i, beta_f, w, gamma, t):

    delta_t, d_delta_t = Delta_Beta(w, beta_i, beta_f, gamma, t)
    ff = Distribution('bose', beta_f, w)
    
    return -1 + (abs(mu)**2)*np.exp(-gamma*t)/ff + delta_t/ff + np.log(ff/delta_t)


def EquidistantInitial(Kinit, beta_eq, w, gamma, T_list, mu_cooling, mu_heating):

    K_heat = lambda beta_c: RelativeEntropy(mu_heating, beta_c, beta_eq, w, gamma, 0) - Kinit
    K_diff = lambda beta_h: ((fc - Distribution('bose', beta_h, w))/ff) + np.log(Distribution('bose', beta_h, w)/fc) - (abs(mu_cooling)**2)/ff

    beta_list = []
    K_heat_list = []
    
    for T in T_list:
        beta = 1/T
        beta_list.append(beta)
        K_heat_list.append(K_heat(beta))

    ## cold
    for i in range(0, np.argmin(K_heat_list), 1):

        if K_heat(beta_list[i]) * K_heat(beta_list[i+1]) < 0:  # houve cruzamento

            beta_1 = brentq(K_heat, beta_list[i], beta_list[i+1])  # raiz exata
            K1 = RelativeEntropy(mu_heating, beta_1, beta_eq, w, gamma, 0)
            break

    ## hot
    fc = Distribution('bose', beta_1, w)
    ff = Distribution('bose', beta_eq, w)

    for i in range(np.argmin(K_heat_list), len(beta_list)-1, 1):

        if K_diff(beta_list[i]) * K_diff(beta_list[i+1]) < 0:  # houve cruzamento

            beta_2 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K2 = RelativeEntropy(mu_cooling, beta_2, beta_eq, w, gamma, 0)
    

    ## plot
    
    Klist = [RelativeEntropy(0, beta_i, beta_eq, w, gamma, 0) for beta_i in beta_list]
    
    plt.plot(T_list, Klist, color='orange')
    plt.scatter([1/beta_eq], [RelativeEntropy(0, beta_eq, beta_eq, w, gamma, 0)], color='orange', label=f'Tw = {1/beta_eq:.3f}')
    plt.scatter([1/beta_1], [K1], color='blue', label=f'Tc = {1/beta_1:.3f}')
    plt.scatter([1/beta_2], [K2], color='red', label=f'Th = {1/beta_2:.3f}')
    plt.hlines(y=Kinit, xmin=min(beta_list), xmax=max(beta_list), color='black', label='Initial Relative Entropy')
    plt.xlabel(r'$\beta$')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()
    
    
    return beta_2, 1/beta_2, K2, beta_1, 1/beta_1, K1





########### MAIN ###########


## Parameters

w = 1
gamma = 0.1
mu_heating = 0
mu_cooling = 1

Kinit = 1

Teq = 2
beta_eq = 1/Teq

T_list = np.arange(0.1, 10, 0.001)

beta_hot, Thot, Khot, beta_cold, Tcold, Kcold = EquidistantInitial(Kinit, beta_eq, w, gamma, T_list, mu_cooling, mu_heating)

tlist = np.arange(0, 100, 0.1)


## Thermal Kinematics

Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating = ThermalKinematics(mu_heating, beta_cold, beta_eq, w, gamma, tlist)
Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling = ThermalKinematics(mu_cooling, beta_hot, beta_eq, w, gamma, tlist)


## Graphs

#### Wigner Fisher Information

plt.plot(tlist, Iw_heating, color='red', label='Heating '+r'$\mu$ = '+f'{mu_heating}')
plt.plot(tlist, Iw_cooling, color='blue', label='Cooling '+r'$\mu$ = '+f'{mu_cooling}')
plt.ylabel('Wigner Fisher Information', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()


#### Velocity

plt.plot(tlist, Vw_heating, color='red', label='Heating '+r'$\mu$ = '+f'{mu_heating}')
plt.plot(tlist, Vw_cooling, color='blue', label='Cooling '+r'$\mu$ = '+f'{mu_cooling}')
plt.ylabel('Statistical Velocity', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()



#### Distance

plt.plot(tlist, Lw_heating, color='red', label='Heating '+r'$\mu$ = '+f'{mu_heating}')
plt.plot(tlist, Lw_cooling, color='blue', label='Cooling '+r'$\mu$ = '+f'{mu_cooling}')
plt.ylabel('Statistical Distance', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()


#### Completion

plt.plot(tlist, completion_heating, color='red', label='Heating '+r'$\mu$ = '+f'{mu_heating}')
plt.plot(tlist, completion_cooling, color='blue', label='Cooling '+r'$\mu$ = '+f'{mu_cooling}')
plt.ylabel('Degree of Completion', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()


#### Relative Entropy

plt.plot(tlist, Kevol_heating, color='red', label='Heating '+r'$\mu$ = '+f'{mu_heating}')
plt.plot(tlist, Kevol_cooling, color='blue', label='Cooling '+r'$\mu$ = '+f'{mu_cooling}')
plt.ylabel('Relative Entropy', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()


#### Entropy Production

plt.plot(tlist, Sprod_heating, color='red', label='Heating '+r'$\mu$ = '+f'{mu_heating}')
plt.plot(tlist, Sprod_cooling, color='blue', label='Cooling '+r'$\mu$ = '+f'{mu_cooling}')
plt.ylabel('Entropy Production', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()

