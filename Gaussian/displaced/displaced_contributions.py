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
    

def Wigner_Fisher_Info(mu, w, beta_i, beta_f, gamma, t):

    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    Iw_passive = 2 * (d_delta_beta/delta_beta)**2
    Iw_ergotropic = 2 * (w**2 + (gamma**2)/4) * (abs(mu)**2)*np.exp(-gamma*t)/delta_beta
    
    return Iw_passive, Iw_ergotropic
    

def EntropyProduction(mu, w, beta_i, beta_f, gamma, t):
    
    ff = Distribution('bose', beta_f, w)
    
    delta_beta, d_delta_beta = Delta_Beta(w, beta_i, beta_f, gamma, t)
    
    Sprod_p = - (d_delta_beta/ff) + (d_delta_beta/delta_beta)
    Sprod_e = (gamma*(abs(mu)**2)*np.exp(-gamma*t))/ff
    
    return Sprod_p, Sprod_e


def ThermalKinematics(mu, beta_i, beta_f, w, gamma, tlist):

    Iw_p = []
    Iw_e = []
    Kevol_p = []
    Kevol_e = []
    Sprod_p = []
    Sprod_e = []

    for t in tlist:

        Iwt_p, Iwt_e = Wigner_Fisher_Info(mu, w, beta_i, beta_f, gamma, t)

        Iw_p.append(Iwt_p)
        Iw_e.append(Iwt_e)

        Kevolt_p, Kevolt_e = RelativeEntropy(mu, beta_i, beta_f, w, gamma, t)
        
        Kevol_p.append(Kevolt_p)
        Kevol_e.append(Kevolt_e)
        
        Sprodt_p, Sprodt_e = EntropyProduction(mu, w, beta_i, beta_f, gamma, t)
        
        Sprod_p.append(Sprodt_p)
        Sprod_e.append(Sprodt_e)
        
    return Iw_p, Iw_e, Kevol_p, Kevol_e, Sprod_p, Sprod_e
    

def RelativeEntropy(mu, beta_i, beta_f, w, gamma, t):

    delta_t, d_delta_t = Delta_Beta(w, beta_i, beta_f, gamma, t)
    ff = Distribution('bose', beta_f, w)
    
    S_p = -1 + delta_t/ff + np.log(ff/delta_t)
    S_e = (abs(mu)**2)*np.exp(-gamma*t)/ff
    
    return S_p, S_e


def EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list):

    K_diff = lambda beta_i: RelativeEntropy(0, beta_i, beta_eq, w, gamma, 0)[0] - Kinit
    
    K_diff_list = []
    for beta in beta_list:
    
        K_diff_list.append(K_diff(beta))

    ## hot
    for i in range(0, np.argmin(K_diff_list), 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_1 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K1 = RelativeEntropy(0, beta_1, beta_eq, w, gamma, 0)[0]

    ## cold
    for i in range(np.argmin(K_diff_list), len(beta_list) - 1, 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_2 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K2 = RelativeEntropy(0, beta_2, beta_eq, w, gamma, 0)[0]
    
    
    ## plot
    
    K = [RelativeEntropy(0, beta_i, beta_eq, w, gamma, 0)[0] for beta_i in beta_list]
    
    plt.plot(beta_list, K, color='orange')
    plt.scatter([beta_eq], [RelativeEntropy(0, beta_eq, beta_eq, w, gamma, 0)[0]], color='orange', label=f'Tw = {1/beta_eq:.3f}')
    plt.scatter([beta_2], [K2], color='blue', label=f'Tc = {1/beta_2:.3f}')
    plt.scatter([beta_1], [K1], color='red', label=f'Th = {1/beta_1:.3f}')
    plt.hlines(y=Kinit, xmin=min(beta_list), xmax=max(beta_list), color='black', label='Initial Relative Entropy')
    plt.xlabel(r'$\beta$')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()
    
    
    return beta_1, 1/beta_1, K1, beta_2, 1/beta_2, K2
        

def WriteOutput(mu, processo, tlist, Iw_p, Iw_e, Kevol_p, Kevol_e, Sprod_p, Sprod_e):        

    f = open(f'./ThermalKinematics/mu{mu}-{processo}-contributions.txt', 'w')
    
    for i in range(len(tlist)):
    
        f.write(f'{tlist[i]} {Iw_p[i]} {Iw_e[i]} {Kevol_p[i]} {Kevol_e[i]} {Sprod_p[i]} {Sprod_e[i]}\n')
    
    f.close()




########### MAIN ###########


## Parameters

w = 1
gamma = 0.1
mu = 0.5

Kinit = 1

Teq = 2
beta_eq = 1/Teq

beta_list = np.arange(0.1, 5, 0.01)

beta_hot, Thot, Khot, beta_cold, Tcold, Kcold = EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list)

tlist = np.arange(0, 100, 0.1)


## Thermal Kinematics

Iw_p_heating, Iw_e_heating, Kevol_p_heating, Kevol_e_heating, Sprod_p_heating, Sprod_e_heating = ThermalKinematics(mu, beta_cold, beta_eq, w, gamma, tlist)
Iw_p_cooling, Iw_e_cooling, Kevol_p_cooling, Kevol_e_cooling, Sprod_p_cooling, Sprod_e_cooling = ThermalKinematics(mu, beta_hot, beta_eq, w, gamma, tlist)


## Write Output Files

WriteOutput(mu, 'heating', tlist, Iw_p_heating, Iw_e_heating, Kevol_p_heating, Kevol_e_heating, Sprod_p_heating, Sprod_e_heating)
WriteOutput(mu, 'cooling', tlist, Iw_p_cooling, Iw_e_cooling, Kevol_p_cooling, Kevol_e_cooling, Sprod_p_cooling, Sprod_e_cooling)










