import numpy as np
import matplotlib.pyplot as plt

def Distribution(beta, dmn):

    f = 1/(np.exp(beta * dmn) - 1)

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
    
def Contributions_Squeezing(gamma, ff, fi, ri, t):
    
    ## delta beta
    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
    
    ## f(beta_t)
    f_beta_t, df_beta_t = F_Beta_t(fi, ff, delta_beta, d_delta_beta, ri, gamma, t)
    
    ## dS_W
    dS = d_Wigner_Entropy(f_beta_t, df_beta_t)        
    
    ## cosh(2r_t)
    cosh_2rt = (delta_beta + 2*fi*np.exp(-gamma*t)*np.sinh(ri)**2)/f_beta_t
    
    ## dr_t
    x = (delta_beta + 2*fi*np.exp(-gamma*t)*np.sinh(ri)**2)/f_beta_t
    dx = ((d_delta_beta - 2*gamma*np.exp(-gamma*t)*fi*np.sinh(ri)**2)*f_beta_t - (delta_beta + 2*np.exp(-gamma*t)*fi*np.sinh(ri)**2)*df_beta_t)/(f_beta_t**2)
    drt = 0.5*dx/np.sqrt((x**2) - 1)
    
    
    ## contribuicao passiva
    pi = dS * (1- (f_beta_t/ff))
    
    ## contribuicao ergotropica
    ergo = -(df_beta_t/ff) * (cosh_2rt - 1) - 2*(f_beta_t/ff) * np.sinh(np.acosh(cosh_2rt)) * drt
    
    return pi, ergo


def Contributions_Displacement(gamma, ff, fi, mu, w, t):
    
    ## delta beta
    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
    
    Sprod_p = d_delta_beta/delta_beta - d_delta_beta/ff
    
    Sprod_e = -(-gamma*(abs(mu)**2)*np.exp(-gamma*t))/(w*ff)
    
    return Sprod_p, Sprod_e


## Parameters

w = 1
gamma = 0.1

Teq = 2
beta_eq = 1/Teq
ff = Distribution(beta_eq, w)

tlist = np.arange(0, 100, 0.1)

muList, beta_hot_list, rList, beta_cold_list = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

symbols = ['-', '--', ':']
labels = ['Total', 'Passive Contribution', 'Ergotropic Contribution']

for i in range(2):
    
    # hot = cooling = displacement || cold = heating = squeezing
    
    fc = Distribution(beta_cold_list[i], w)
    fh = Distribution(beta_hot_list[i], w)
    
    r = rList[i]
    mu = muList[i]

    # ergotropia
    ergo_cooling_list = []
    ergo_heating_list = []
    
    # producao entropia passivo
    pi_cooling_list = []
    pi_heating_list = []
    
    for t in tlist:
        
        pi_cooling, ergo_cooling = Contributions_Displacement(gamma, ff, fh, mu, w, t)
        pi_heating, ergo_heating = Contributions_Squeezing(gamma, ff, fc, r, t)

        pi_cooling_list.append(pi_cooling)
        pi_heating_list.append(pi_heating)
        
        ergo_cooling_list.append(ergo_cooling)
        ergo_heating_list.append(ergo_heating)
        
    
    Sprod_cooling = np.loadtxt(f'./ThermalKinematics/d{mu}-cooling.txt', unpack=True, usecols=(6))
    Sprod_heating = np.loadtxt(f'./ThermalKinematics/s{r}-heating.txt', unpack=True, usecols=(6))
    
    
    fig = plt.figure(figsize=(12,6))

    fig_heating = plt.subplot(1, 2, 1)
    
    line_ht = plt.plot(tlist, Sprod_heating, color='red', linestyle=symbols[0], linewidth=2)
    line_hp = plt.plot(tlist, pi_heating_list, color='red', linestyle=symbols[1], linewidth=2)
    line_he = plt.plot(tlist, ergo_heating_list, color='red', linestyle=symbols[2], linewidth=2)
    
    plt.xscale('log')
    plt.xlim(left=0.1)
    plt.title(r'$r$ = '+f'{rList[i]:.2f}', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.ylabel('Entropy Production Rate', fontsize=12)
    handles_red = [line_ht[0], line_hp[0], line_he[0]]
    
    
    fig_cooling = plt.subplot(1, 2, 2)

    line_ct = plt.plot(tlist, Sprod_cooling, color='blue', linestyle=symbols[0], linewidth=2)
    line_cp = plt.plot(tlist, pi_cooling_list, color='blue', linestyle=symbols[1], linewidth=2)
    line_ce = plt.plot(tlist, ergo_cooling_list, color='blue', linestyle=symbols[2], linewidth=2)

    plt.xscale('log')
    plt.xlim(left=0.1)
    plt.title(r'$\mu$ = '+f'{muList[i]:.2f}', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    handles_blue = [line_ct[0], line_cp[0], line_ce[0]]


    # Legenda Heating (Vermelha) - Superior
    leg_h = fig.legend(handles_red, labels, 
                       loc='lower center', 
                       ncol=3, 
                       title="Heating with Squeezing", 
                       title_fontproperties={'weight':'bold', 'size':12},
                       fontsize=12,
                       bbox_to_anchor=(0.5, 0.1), 
                       frameon=False)

    # Legenda Cooling (Azul) - Inferior
    leg_c = fig.legend(handles_blue, labels, 
                       loc='lower center', 
                       ncol=3, 
                       title="Cooling with Displacement", 
                       title_fontproperties={'weight':'bold', 'size':12},
                       fontsize=12,
                       bbox_to_anchor=(0.5, 0.01), 
                       frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.35)
    plt.show()


















