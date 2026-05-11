import numpy as np
import matplotlib.pyplot as plt

def Distribution(beta, dmn):

    f = 1/(np.exp(beta * dmn) - 1)

    return f


## Parameters

w = 1
gamma = 0.1

Teq = 2
beta_eq = 1/Teq

tlist = np.arange(0, 50, 0.1)

r_list, beta_c_list, beta_h_list = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

f = plt.figure(figsize=(5,10))
symbols = ['-', '--']

for i, r in enumerate(r_list):

    fc = Distribution(beta_c_list[i], w)
    fh = Distribution(beta_h_list[i], w)
    ff = Distribution(beta_eq, w)
    
    # ergotropia
    ergo_c = []
    ergo_h = []
    
    # passivo
    pi_c = []
    pi_h = []
    
    # razao ergotropia / passivo
    Rt_c = []
    Rt_h = []
    
    print(fh - fc + ff*np.log(fc/fh))
    
    for t in tlist:
    
        delta_beta_c = (fc - ff) * np.exp(-gamma*t) + ff
        delta_beta_h = (fh - ff) * np.exp(-gamma*t) + ff
        
        f_beta_t_c = np.sqrt( delta_beta_c**2 + 4*fc*ff*np.exp(-2*gamma*t)*(np.exp(gamma*t) - 1)*np.sinh(r)**2 )
        f_beta_t_h = np.sqrt( delta_beta_h**2 + 4*fh*ff*np.exp(-2*gamma*t)*(np.exp(gamma*t) - 1)*np.sinh(r)**2 )
        
        cosh_2rt_c = (delta_beta_c + 2*fc*np.exp(-gamma*t)*np.sinh(r)**2)/f_beta_t_c
        cosh_2rt_h = (delta_beta_h + 2*fh*np.exp(-gamma*t)*np.sinh(r)**2)/f_beta_t_h

        pi_c.append(w*f_beta_t_c)
        pi_h.append(w*f_beta_t_h)
        
        ergo_c.append(w*f_beta_t_c*(cosh_2rt_c - 1))
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























