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


def EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list):

    K_diff = lambda beta_i: RelativeEntropy(0, beta_i, beta_eq, w, gamma, 0) - Kinit
    print(RelativeEntropy(0, beta_eq, beta_eq, w, gamma, 0))
    K_diff_list = []
    for beta in beta_list:
    
        K_diff_list.append(K_diff(beta))

    ## hot
    for i in range(0, np.argmin(K_diff_list), 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_1 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K1 = RelativeEntropy(0, beta_1, beta_eq, w, gamma, 0)

    ## cold
    for i in range(np.argmin(K_diff_list), len(beta_list) - 1, 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_2 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            K2 = RelativeEntropy(0, beta_2, beta_eq, w, gamma, 0)
    
    
    ## plot
    
    K = [RelativeEntropy(0, beta_i, beta_eq, w, gamma, 0) for beta_i in beta_list]
    
    plt.plot(beta_list, K, color='orange')
    plt.scatter([beta_eq], [RelativeEntropy(0, beta_eq, beta_eq, w, gamma, 0)], color='orange', label=f'Tw = {1/beta_eq:.3f}')
    plt.scatter([beta_2], [K2], color='blue', label=f'Tc = {1/beta_2:.3f}')
    plt.scatter([beta_1], [K1], color='red', label=f'Th = {1/beta_1:.3f}')
    plt.hlines(y=Kinit, xmin=min(beta_list), xmax=max(beta_list), color='black', label='Initial Relative Entropy')
    plt.xlabel(r'$\beta$')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()
    
    
    return beta_1, 1/beta_1, K1, beta_2, 1/beta_2, K2
        

def WriteOutput(mu, processo, tlist, Iw, Vw, Lw, completion, Kevol, Sprod):        

    f = open(f'./ThermalKinematics/mu{mu}-{processo}.txt', 'w')
    
    for i in range(len(tlist)):
    
        f.write(f'{tlist[i]} {Iw[i]} {Vw[i]} {Lw[i]} {completion[i]} {Kevol[i]} {Sprod[i]}\n')
    
    f.close()




########### MAIN ###########


## Parameters

w = 1
gamma = 0.1
mu = 2

Kinit = 1

Teq = 2
beta_eq = 1/Teq

beta_list = np.arange(0.1, 5, 0.01)

beta_hot, Thot, Khot, beta_cold, Tcold, Kcold = EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list)

tlist = np.arange(0, 100, 0.1)


## Thermal Kinematics

Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating = ThermalKinematics(mu, beta_cold, beta_eq, w, gamma, tlist)
Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling = ThermalKinematics(mu, beta_hot, beta_eq, w, gamma, tlist)


## Write Output Files

WriteOutput(mu, 'heating', tlist, Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating)
WriteOutput(mu, 'cooling', tlist, Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling)










