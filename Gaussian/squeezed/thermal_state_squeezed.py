import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def Distribution(dist, beta_R, dmn):

    if dist == 'bose':
        f = 1/(np.exp(beta_R * dmn) - 1)

    elif dist == 'fermi':
        f = 1/(np.exp(beta_R * dmn) + 1)
        
    return f


def Delta_Beta(fi, ff, gamma, t):

    delta_beta = (fi - ff) * np.exp(-gamma*t) + ff
    
    d_delta_beta = -gamma * (fi - ff) * np.exp(-gamma*t)
    
    return delta_beta, d_delta_beta


def F_Beta_t(fi, ff, delta_beta, d_delta_beta, r, gamma, t):

    f = np.sqrt( delta_beta**2 + 4*fi*ff*np.exp(-2*gamma*t)*(np.exp(gamma*t) - 1)*np.sinh(r)**2 )
    
    df = (0.5/f) * ( 2*delta_beta*d_delta_beta + 4*fi*ff*gamma*np.exp(-gamma*t)*(2*np.exp(-gamma*t) - 1)*np.sinh(r)**2)
    
    return f, df
    

def d_Wigner_Entropy(fb, dfb):
    
    dSw = dfb/fb
    
    return dSw


def d_Rt(com_derivada, fi, fb, dfb, delta_beta, d_delta_beta, r, gamma, t):

    x = (delta_beta + 2*fi*np.exp(-gamma*t)*np.sinh(r)**2)/fb
    
    rt = 0.5*np.arccosh(x)
    
    if com_derivada:
        dx = ((d_delta_beta - 2*gamma*np.exp(-gamma*t)*fi*np.sinh(r)**2)*fb - (delta_beta + 2*np.exp(-gamma*t)*fi*np.sinh(r)**2)*dfb)/(fb**2)
        
        drt = -0.5*dx/np.sqrt(1 - x**2)
    
        return rt, drt
    
    else:
        return rt
    

def EntropyProduction(ff, fb, dfb, rt, drt, dSw):
    
    dE = w*(dfb*np.cosh(2*rt) + fb*np.sinh(2*rt)*2*drt)
    
    Sprod = dSw - dE/(w*ff)
    
    return Sprod


def RelativeEntropy(fi, ff, r, gamma, t):

    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
    
    fb, dfb = F_Beta_t(fi, ff, delta_beta, d_delta_beta, r, gamma, t)
    
    rt = d_Rt(False, fi, fb, dfb, delta_beta, d_delta_beta, r, gamma, t)
    
    return -1 + (fb/ff)*np.cosh(2*rt) + np.log(ff/fb)
    

def Wigner_Fisher_Info(rt, drt, w, dSw):
    
    dphi = -2*w
    
    return 2*(dSw**2) + 8*drt**2 + 2*dphi*np.sinh(2*rt)**2
    

def Velocity(Iw):

    return np.sqrt(Iw)
    

def Distance(vList, dt):
    
    L = np.cumsum(vList) * dt
    
    return L


def ThermalKinematics(r, beta_i, beta_f, w, gamma, tlist):

    Iw = []
    Vw = []
    Kevol = []
    Sprod = []

    fi = Distribution('bose', beta_i, w)
    ff = Distribution('bose', beta_f, w)

    for t in tlist:
        
        delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
        
        fb, dfb = F_Beta_t(fi, ff, delta_beta, d_delta_beta, r, gamma, t)

        rt, drt = d_Rt(True, fi, fb, dfb, delta_beta, d_delta_beta, r, gamma, t)
        
        dSw = d_Wigner_Entropy(fb, dfb)
    
    
        Iwt = Wigner_Fisher_Info(rt, drt, w, dSw)

        Vwt = Velocity(Iwt)

        Iw.append(Iwt)
        Vw.append(Vwt)

        Kevol.append(RelativeEntropy(fi, ff, r, gamma, t))
        
        Sprod.append(EntropyProduction(ff, fb, dfb, rt, drt, dSw))
    
    Lw = Distance(Vw, tlist[1]-tlist[0]) 

    completion = []

    for i in range(len(Lw)):

        completion.append(Lw[i]/Lw[-1])
        
    return Iw, Vw, Lw, completion, Kevol, Sprod
    

def EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list):

    ff = Distribution('bose', beta_eq, w)

    K_diff = lambda fi: RelativeEntropy(fi, ff, 0, gamma, 0) - Kinit

    K_diff_list = []
    fi_list = []
    
    for beta_i in beta_list:
        fi = Distribution('bose', beta_i, w)
        fi_list.append(fi)
        K_diff_list.append(K_diff(fi))

    ## hot
    for i in range(0, np.argmin(K_diff_list), 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_1 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            f1 = Distribution('bose', beta_1, w)
            K1 = RelativeEntropy(f1, ff, 0, gamma, 0)

    ## cold
    for i in range(np.argmin(K_diff_list), len(beta_list) - 1, 1):

        if K_diff_list[i] * K_diff_list[i+1] < 0:  # houve cruzamento

            beta_2 = brentq(K_diff, beta_list[i], beta_list[i+1])  # raiz exata
            f2 = Distribution('bose', beta_2, w)
            K2 = RelativeEntropy(f2, ff, 0, gamma, 0)
    
    
    ## plot
    
    K = [RelativeEntropy(fi, ff, 0, gamma, 0) for fi in fi_list]
    
    plt.plot(beta_list, K, color='orange')
    plt.scatter([beta_eq], [RelativeEntropy(ff, ff, 0, gamma, 0)], color='orange', label=f'Tw = {1/beta_eq:.3f}')
    #plt.scatter([beta_2], [K2], color='blue', label=f'Tc = {1/beta_2:.3f}')
    #plt.scatter([beta_1], [K1], color='red', label=f'Th = {1/beta_1:.3f}')
    plt.hlines(y=Kinit, xmin=min(beta_list), xmax=max(beta_list), color='black', label='Initial Relative Entropy')
    plt.xlabel(r'$\beta$')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()
    
    
    return beta_1, 1/beta_1, K1, beta_2, 1/beta_2, K2
        

def WriteOutput(r, processo, tlist, Iw, Vw, Lw, completion, Kevol, Sprod):        

    f = open(f'./ThermalKinematics/r{r}-{processo}.txt', 'w')
    
    for i in range(len(tlist)):
    
        f.write(f'{tlist[i]} {Iw[i]} {Vw[i]} {Lw[i]} {completion[i]} {Kevol[i]} {Sprod[i]}\n')
    
    f.close()




########### MAIN ###########


## Parameters

w = 1
gamma = 0.1
r = 1

Kinit = 1

Teq = 2
beta_eq = 1/Teq

beta_list = np.arange(0.1, 10, 0.01)

beta_hot, Thot, Khot, beta_cold, Tcold, Kcold = EquidistantInitial(Kinit, beta_eq, w, gamma, beta_list)

tlist = np.arange(0, 100, 0.1)


## Thermal Kinematics

Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating = ThermalKinematics(r, beta_cold, beta_eq, w, gamma, tlist)
Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling = ThermalKinematics(r, beta_hot, beta_eq, w, gamma, tlist)


## Write Output Files

WriteOutput(r, 'heating', tlist, Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating)
WriteOutput(r, 'cooling', tlist, Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling)










