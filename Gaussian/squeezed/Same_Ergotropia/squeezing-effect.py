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
    
def Contributions(dS, f_beta_t, ff, fi, ri):
    
    cosh_2rt = (delta_beta + 2*fi*np.exp(-gamma*t)*np.sinh(ri)**2)/f_beta_t
    
    pi = dS * (1- (f_beta_t/ff))
    ergo = w*f_beta_t*(cosh_2rt - 1)
    Rt = ergo/pi
    
    return pi, ergo, Rt


## Parameters

w = 1
gamma = 0.1

Teq = 2
beta_eq = 1/Teq
ff = Distribution(beta_eq, w)

tlist = np.arange(0, 100, 0.1)

rlist_cooling, beta_list_cooling, rlist_heating, beta_list_heating = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

f = plt.figure(figsize=(5,9))
symbols = ['-', '--']

for i in range(2)

    fc = Distribution(beta_list_heating[i], w)
    fh = Distribution(beta_list_cooling[i], w)
    
    # ergotropia
    ergo_c = []
    ergo_h = []
    
    # passivo
    pi_c = []
    pi_h = []
    
    # razao ergotropia / passivo
    Rt_c = []
    Rt_h = []
    
    for t in tlist:
    
        delta_beta_c, d_delta_beta_c = Delta_Beta(fc, ff, gamma, t)
        delta_beta_h, d_delta_beta_h = Delta_Beta(fh, ff, gamma, t)
        
        f_beta_t_c, df_beta_t_c = F_Beta_t(fc, ff, delta_beta_c, d_delta_beta_c, rlist_heating[i], gamma, t)
        f_beta_t_h, df_beta_t_h = F_Beta_t(fh, ff, delta_beta_h, d_delta_beta_h, rlist_cooling[i], gamma, t)
        
        dSc = d_Wigner_Entropy(f_beta_t_c, df_beta_t_c)
        dSh = d_Wigner_Entropy(f_beta_t_h, df_beta_t_h)

        pi_c.append()
        pi_h.append(dSh * (1- (f_beta_t_h/ff)))
        
        ergo_c.append()
        ergo_h.append(w*f_beta_t_h*(cosh_2rt_h - 1))
        
        Rt_c.append(cosh_2rt_c - 1)
        Rt_h.append(cosh_2rt_h - 1)
   
    plt.subplot(311)
    plt.plot(tlist, pi_c, color='red', linestyle=symbols[i], label=f'Heating - r = {r}')
    plt.plot(tlist, pi_h, color='blue', linestyle=symbols[i], label=f'Cooling - r = {r}')
    
    plt.subplot(312)
    plt.plot(tlist, ergo_c, color='red', linestyle=symbols[i], label=f'Heating - r = {r}')
    plt.plot(tlist, ergo_h, color='blue', linestyle=symbols[i], label=f'Cooling - r = {r}')
   
    plt.subplot(313)
    plt.plot(tlist, Rt_c, color='red', linestyle=symbols[i], label=f'Heating - r = {r}')
    plt.plot(tlist, Rt_h, color='blue', linestyle=symbols[i], label=f'Cooling - r = {r}')


plt.subplot(311)
plt.xlabel('Time')
plt.ylabel(r'$E_{\pi}$')
plt.legend()
plt.title('Passive State Energy')

plt.subplot(312)
plt.xlabel('Time')
plt.ylabel(r'$\mathcal{E}$')
plt.legend()
plt.title('Ergotropy')

plt.subplot(313)
plt.xlabel('Time')
plt.ylabel(r'$R = \frac{\mathcal{E}}{E_{\pi}}$')
plt.legend()
plt.title('Ratio between Ergotropy and Passive State Energy')


plt.tight_layout()
plt.show()























