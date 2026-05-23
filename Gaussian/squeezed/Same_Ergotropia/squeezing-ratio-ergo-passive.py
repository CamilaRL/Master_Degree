import numpy as np
import matplotlib.pyplot as plt

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
    
    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
        
    f_beta_t, df_beta_t = F_Beta_t(fi, ff, delta_beta, d_delta_beta, ri, gamma, t)
        
    dS = d_Wigner_Entropy(f_beta_t, df_beta_t)        
    
    cosh_2rt = (delta_beta + 2*fi*np.exp(-gamma*t)*np.sinh(ri)**2)/f_beta_t
    
    E_pi = w*f_beta_t
    ergo = E_pi*(cosh_2rt - 1)
    Rt = ergo/E_pi
    
    return E_pi, ergo, Rt


## Parameters

w = 1
gamma = 0.1

Teq = 4
beta_eq = 1/Teq
ff = Distribution(beta_eq, w)

tlist = np.arange(0, 40, 0.1)

rlist_hot, beta_list_hot, rlist_cool, beta_list_cool = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

f = plt.figure(figsize=(12,5))
symbols = ['-', '--']

for i in range(2):

    fc = Distribution(beta_list_cool[i], w)
    fh = Distribution(beta_list_hot[i], w)
    
    rc = rlist_cool[i]
    rh = rlist_hot[i]

    # ergotropia
    ergo_c_list = []
    ergo_h_list = []
    
    # producao entropia passivo
    Epi_c_list = []
    Epi_h_list = []
    
    # razao ergotropia / energia passivo
    Rt_c_list = []
    Rt_h_list = []
    
    for t in tlist:
        
        Epi_c, ergo_c, Rt_c = Contributions(gamma, ff, fc, rc, t)
        Epi_h, ergo_h, Rt_h = Contributions(gamma, ff, fh, rh, t)

        Epi_c_list.append(Epi_c)
        Epi_h_list.append(Epi_h)
        
        ergo_c_list.append(ergo_c)
        ergo_h_list.append(ergo_h)
        
        Rt_c_list.append(Rt_c)
        Rt_h_list.append(Rt_h)
   
    plt.subplot(131)
    plt.plot(tlist, Epi_c_list, color='red', linestyle=symbols[i], linewidth=2, label=f'Heating - r = {rc:.2f}')
    plt.plot(tlist, Epi_h_list, color='blue', linestyle=symbols[i], linewidth=2, label=f'Cooling - r = {rh:.2f}')
    
    plt.subplot(132)
    plt.plot(tlist, ergo_c_list, color='red', linestyle=symbols[i], linewidth=2, label=f'Heating - r = {rc:.2f}')
    plt.plot(tlist, ergo_h_list, color='blue', linestyle=symbols[i], linewidth=2, label=f'Cooling - r = {rh:.2f}')
   
    fig3 = plt.subplot(133)
    plt.plot(tlist, Rt_c_list, color='red', linestyle=symbols[i], linewidth=2, label=f'Heating - r = {rc:.2f}')
    plt.plot(tlist, Rt_h_list, color='blue', linestyle=symbols[i], linewidth=2, label=f'Cooling - r = {rh:.2f}')


plt.subplot(131)
plt.xlabel('Time', fontsize=12)
plt.ylabel(r'$E_{\pi}$', fontsize=12)
plt.title('Passive State Energy', fontsize=14)

plt.subplot(132)
plt.xlabel('Time', fontsize=12)
plt.ylabel(r'$\mathcal{E}$', fontsize=12)
plt.title('Ergotropy', fontsize=14)

plt.subplot(133)
plt.xlabel('Time', fontsize=12)
plt.ylabel(r'$R = \frac{\mathcal{E}}{E_{\pi}}$', fontsize=12)
plt.title(f'Ratio between Ergotropy and \nPassive State Energy', fontsize=14)

handles, labels = fig3.get_legend_handles_labels()

plt.tight_layout()
plt.legend(handles, labels, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(-0.75, -0.35), fontsize=12)
plt.subplots_adjust(bottom=0.25)
plt.show()























