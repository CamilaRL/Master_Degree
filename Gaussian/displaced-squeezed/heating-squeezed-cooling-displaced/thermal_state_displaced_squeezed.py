import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def Distribution(beta_R, dmn):

    n = 1/(np.exp(beta_R * dmn) - 1)
    
    f = n + 0.5
        
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

        drt = 0.5*dx/np.sqrt((x**2) - 1)
            
        return rt, drt
    
    else:
        return rt


def EntropyProduction_Squeezing(ff, fb, dfb, rt, drt, dSw):
    
    dE = w*(dfb*np.cosh(2*rt) + fb*np.sinh(2*rt)*2*drt)
    
    Sprod = dSw - dE/(w*ff)
    
    return Sprod


def EntropyProduction_Displacement(mu, delta_beta, d_delta_beta, ff, gamma, t):
    
    Sprod = (gamma*(abs(mu)**2)*np.exp(-gamma*t) - d_delta_beta)/ff + (d_delta_beta/delta_beta)
    
    return Sprod


def RelativeEntropy_Squeezing(r, fi, ff, gamma, t):

    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
    
    fb, dfb = F_Beta_t(fi, ff, delta_beta, d_delta_beta, r, gamma, t)
    
    rt = d_Rt(False, fi, fb, dfb, delta_beta, d_delta_beta, r, gamma, t)
    
    return -1 + (fb/ff)*np.cosh(2*rt) + np.log(ff/fb)
    

def RelativeEntropy_Displacement(mu, delta_t, ff, gamma, t):

    return -1 + (abs(mu)**2)*np.exp(-gamma*t)/ff + delta_t/ff + np.log(ff/delta_t)


def Wigner_Fisher_Info_Squeezing(rt, drt, w, dSw):
    
    dphi = -2*w

    return 2*(dSw**2) + 8*drt**2 + 2*(dphi*np.sinh(2*rt))**2


def Wigner_Fisher_Info_Displacement(mu, w, delta_beta, d_delta_beta, gamma, t):
    
    return 2 * (d_delta_beta/delta_beta)**2 + 2 * (w**2 + (gamma**2)/4) * (abs(mu)**2)*np.exp(-gamma*t)/delta_beta
    

def Velocity(Iw):

    return np.sqrt(Iw/2)
    

def Distance(vList, dt):
    
    L = np.cumsum(vList) * dt
    
    return L


def ThermalKinematics(recurso, param, beta_i, beta_f, w, gamma, tlist):

    Iw = []
    Vw = []
    Kevol = []
    Sprod = []

    fi = Distribution(beta_i, w)
    ff = Distribution(beta_f, w)

    if recurso == 's':

        for t in tlist:
        
            delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)

            fb, dfb = F_Beta_t(fi, ff, delta_beta, d_delta_beta, param, gamma, t)

            rt, drt = d_Rt(True, fi, fb, dfb, delta_beta, d_delta_beta, param, gamma, t)

            dSw = d_Wigner_Entropy(fb, dfb)
        
            Kevol.append(RelativeEntropy_Squeezing(r, fi, ff, gamma, t))
            
            Sprod.append(EntropyProduction_Squeezing(ff, fb, dfb, rt, drt, dSw))
    
            Iwt = Wigner_Fisher_Info_Squeezing(rt, drt, w, dSw)
            
            Vwt = Velocity(Iwt)

            Iw.append(Iwt)
            Vw.append(Vwt)


    elif recurso == 'd':
        
        for t in tlist:
        
            delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
            
            Kevol.append(RelativeEntropy_Displacement(param, delta_beta, ff, gamma, t))
            
            Sprod.append(EntropyProduction_Displacement(param, delta_beta, d_delta_beta, ff, gamma, t))
            
            Iwt = Wigner_Fisher_Info_Displacement(param, w, delta_beta, d_delta_beta, gamma, t)

            Vwt = Velocity(Iwt)

            Iw.append(Iwt)
            Vw.append(Vwt)

    
    Lw = Distance(Vw, tlist[1]-tlist[0]) 

    completion = []

    for i in range(len(Lw)):

        completion.append(Lw[i]/Lw[-1])
        
    return Iw, Vw, Lw, completion, Kevol, Sprod
    

def EquidistantInitial(Kinit, beta_eq, mu, w, gamma, beta_list):

    ff = Distribution(beta_eq, w)
    
    Kd_func = lambda fi: RelativeEntropy_Displacement(mu, fi, ff, gamma, 0) - Kinit
    
    Kd_list = []
    fi_list = []
    
    for beta_i in beta_list:
        fi = Distribution(beta_i, w)
        fi_list.append(fi)
        Kd_list.append(Kd_func(fi))
    
    ## displacement = cooling
    for i in range(0, np.argmin(Kd_list), 1):
        
        if Kd_list[i] * Kd_list[i+1] < 0:  # houve cruzamento

            fd = brentq(Kd_func, fi_list[i], fi_list[i+1])  # raiz exata
            Kd = RelativeEntropy_Displacement(mu, fd, ff, gamma, 0)
            beta_d = np.log((1/(fd - 0.5)) + 1)/w
    
    ## squeezing = heating
    Ks_func = lambda fi: fd - fi + ff*np.log(fi/fd)
    
    Ks_list = []
    for fi in fi_list:
        Ks_list.append(Ks_func(fi))
    
    for i in range(np.argmin(Kd_list), len(beta_list) - 1, 1):

        if Ks_list[i] * Ks_list[i+1] < 0:  # houve cruzamento

            fs = brentq(Ks_func, fi_list[i], fi_list[i+1])  # raiz exata
            
            r = 0.5 * np.arccosh((abs(mu)**2)/fs +1)
            
            Ks = RelativeEntropy_Squeezing(r, fs, ff, gamma, 0)
            beta_s = np.log((1/(fs - 0.5)) + 1)/w
    
    
    ## plot
    
    plt.hlines(y=Kinit, xmin=min(beta_list), xmax=max(beta_list), color='grey', label='Initial Relative Entropy')
    plt.scatter([beta_eq], [RelativeEntropy_Displacement(0, ff, ff, gamma, 0)], color='black', label=f'Tw = {1/beta_eq:.3f}')

    Kd_curve = [RelativeEntropy_Displacement(mu, fi, ff, gamma, 0) for fi in fi_list]
    plt.plot(beta_list, Kd_curve, color='red')
    
    Ks_curve = [RelativeEntropy_Squeezing(r, fi, ff, gamma, 0) for fi in fi_list]
    plt.plot(beta_list, Ks_curve, color='blue')
    
    plt.scatter([beta_d], [Kd], color='red', label=f'Td = {1/beta_d:.3f}')
    plt.scatter([beta_s], [Ks], color='blue', label=f'Ts = {1/beta_s:.3f}')
    plt.xlabel(r'$\beta$')
    plt.ylabel('Relative Entropy')
    plt.legend()
    plt.show()
    
    
    return beta_d, 1/beta_d, Kd, beta_s, 1/beta_s, Ks, r


def WriteOutput(filename, tlist, Iw, Vw, Lw, completion, Kevol, Sprod):        

    f = open(f'./ThermalKinematics/{filename}', 'w')
    
    for i in range(len(tlist)):
    
        f.write(f'{tlist[i]} {Iw[i]} {Vw[i]} {Lw[i]} {completion[i]} {Kevol[i]} {Sprod[i]}\n')
    
    f.close()




########### MAIN ###########


## Parameters

w = 1
gamma = 0.1
mu = 1.0

Kinit = 1

Teq = 4
beta_eq = 1/Teq

beta_list = np.arange(0.01, 6, 0.001)

beta_hot, Thot, Khot, beta_cold, Tcold, Kcold, r = EquidistantInitial(Kinit, beta_eq, mu, w, gamma, beta_list)
# hot = cooling = displacement || cold = heating = squeezing

file_init = open('./ThermalKinematics/initial_temperatures.txt', 'a')

file_init.write(f'{mu} {beta_hot} {r} {beta_cold}\n')

file_init.close()


tlist = np.arange(0, 100, 0.1)


## Thermal Kinematics

Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating = ThermalKinematics('s', r, beta_cold, beta_eq, w, gamma, tlist)
Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling = ThermalKinematics('d', mu, beta_hot, beta_eq, w, gamma, tlist)


## Write Output Files

WriteOutput(f's{r}-heating.txt', tlist, Iw_heating, Vw_heating, Lw_heating, completion_heating, Kevol_heating, Sprod_heating)
WriteOutput(f'd{mu}-cooling.txt', tlist, Iw_cooling, Vw_cooling, Lw_cooling, completion_cooling, Kevol_cooling, Sprod_cooling)







