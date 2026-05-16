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


def EntropyProduction(w, beta_i, beta_f, gamma, t):
    
    ff = Distribution('bose', beta_f, w)
    
    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    #Sprod = - (d_delta_beta/ff) + (d_delta_beta/delta_beta)
    Sprod = gamma*((delta_beta/ff) + (ff/delta_beta) - 2)
    
    return Sprod
    

def EntropyFluxRate(w, beta_i, beta_f, gamma, t):
    
    ff = Distribution('bose', beta_f, w)
    
    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    return gamma*((delta_beta/ff) - 1)


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
left_heating = []
left_cooling = []

right_heating = []
right_cooling = []

for t in tlist:

    Sprod_heating = EntropyProduction(w, beta_cold, beta_eq, gamma, t)
    Sprod_cooling = EntropyProduction(w, beta_hot, beta_eq, gamma, t)
    
    flux_heating = EntropyFluxRate(w, beta_cold, beta_eq, gamma, t)
    flux_cooling = EntropyFluxRate(w, beta_hot, beta_eq, gamma, t)
    
    left_heating.append((1/Sprod_heating) - (1/Sprod_cooling))
    right_heating.append(2/flux_heating)
    
    left_cooling.append((1/Sprod_cooling) - (1/Sprod_heating))
    right_cooling.append(2/flux_cooling)





plt.plot(tlist[:50], left_heating[:50], color='red')
plt.plot(tlist[:50], right_heating[:50], color='black', label=r'$\frac{2}{\phi}$')
plt.legend()
plt.show()


plt.plot(tlist, left_cooling, color='blue')
plt.plot(tlist, right_cooling, color='black', label=r'$\frac{2}{\phi}$')
plt.legend()
plt.show()