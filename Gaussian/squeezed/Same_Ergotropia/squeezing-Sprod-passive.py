import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def Distribution(beta, dmn):

    n = 1/(np.exp(beta * dmn) - 1)
    
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
    
def Contributions(gamma, ff, fi, ri, t):
    
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


## Parameters

w = 1
gamma = 0.1

Teq = 4
beta_eq = 1/Teq
ff = Distribution(beta_eq, w)

tlist = np.arange(0, 40, 0.1)

rlist_hot, beta_list_hot, rlist_cool, beta_list_cool = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

symbols = ['-', '--', ':']
labels = ['Total', 'Passive Contribution', 'Ergotropic Contribution']

fig = plt.figure(figsize=(12,5))

for i in range(0,2,1):

    fc = Distribution(beta_list_cool[i], w)
    fh = Distribution(beta_list_hot[i], w)
    
    rc = rlist_cool[i]
    rh = rlist_hot[i]

    # ergotropia
    ergo_c_list = []
    ergo_h_list = []
    
    # producao entropia passivo
    pi_c_list = []
    pi_h_list = []
    
    for t in tlist:
        
        pi_c, ergo_c = Contributions(gamma, ff, fc, rc, t)
        pi_h, ergo_h = Contributions(gamma, ff, fh, rh, t)
        
        '''if pi_c < 0:
            print(f'Heating non-Markovian at time {t} - {rc}')
        if pi_h < 0:
            print(f'Cooling non-Markovian at time {t} - {rh}')'''

        pi_c_list.append(pi_c)
        pi_h_list.append(pi_h)
        
        ergo_c_list.append(ergo_c)
        ergo_h_list.append(ergo_h)
    
    
    print(pi_c_list[0]-pi_h_list[0])
    # Legenda dentro
    
    n = 121 + i

    fig_heating = plt.subplot(n)
    
    plt.plot(tlist, pi_c_list, color='red', linewidth=2)
    plt.plot(tlist, pi_h_list, color='blue', linewidth=2)
    
    plt.xscale('log')
    plt.ylim(top=0.6)
    plt.xlim(left=0.1)
    plt.title(r'$r_h$ = '+f'{rlist_hot[i]:.2f} -'+r'$r_c$ = '+f'{rlist_cool[i]:.2f}', fontsize=16)
    plt.xlabel('Time', fontsize=14)
    plt.xticks(fontsize=14)

    if i == 0:
        plt.yticks(fontsize=14)
        plt.ylabel('Entropy Production Rate', fontsize=14)
    else:
        plt.yticks([])
    
plt.tight_layout()
plt.show()
    





















