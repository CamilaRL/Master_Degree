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
    
def Passive_Squeezing(gamma, ff, fi, ri, t):
    
    ## delta beta
    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
    
    ## f(beta_t)
    f_beta_t, df_beta_t = F_Beta_t(fi, ff, delta_beta, d_delta_beta, ri, gamma, t)
    
    ## dS_W
    dSt = d_Wigner_Entropy(f_beta_t, df_beta_t)
    
    Ipt = 2*(dSt**2)
    
    Vpt = np.sqrt(Ipt)
    
    return Ipt, Vpt


def Passive_Displacement(gamma, ff, fi, mu, w, t):
    
    ## delta beta
    delta_beta, d_delta_beta = Delta_Beta(fi, ff, gamma, t)
    
    dSt = d_delta_beta/delta_beta
    
    Ipt = 2*(dSt**2)
    
    Vpt = np.sqrt(Ipt)
    
    return Ipt, Vpt


def Distance(vList, dt):
    
    L = np.cumsum(vList) * dt
    
    return L

## Parameters

w = 1
gamma = 0.1

Teq = 4

beta_eq = 1/Teq
ff = Distribution(beta_eq, w)

tlist = np.arange(0, 110, 0.1)

rList, beta_hot_list, muList, beta_cold_list = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

symbols = ['-', '--']

for i in range(2):
    
    # hot = cooling = squeezing || cold = heating = displacement
    
    fhot = Distribution(beta_hot_list[i], w)
    fcold = Distribution(beta_cold_list[i], w)
    
    r = rList[i]
    mu = muList[i]

    Ip_squeezing = []
    Vp_squeezing = []
    
    Ip_displace = []
    Vp_displace = []
    
    for t in tlist:
        
        Ipt_squeezing, Vpt_squeezing = Passive_Squeezing(gamma, ff, fhot, r, t)
        Ipt_displace, Vpt_displace = Passive_Displacement(gamma, ff, fcold, mu, w, t)

        Ip_squeezing.append(Ipt_squeezing)
        Vp_squeezing.append(Vpt_squeezing)
        
        Ip_displace.append(Ipt_displace)
        Vp_displace.append(Vpt_displace)
    
    
    Lp_squeezing = Distance(Vp_squeezing, tlist[1]-tlist[0])    
    Lp_displace = Distance(Vp_displace, tlist[1]-tlist[0])
    
    completion_p_squeezing = []
    completion_p_displace = []

    for t in range(len(Lp_displace)):

        completion_p_squeezing.append(Lp_squeezing[t]/Lp_squeezing[-1])
        completion_p_displace.append(Lp_displace[t]/Lp_displace[-1])
    
    
    plt.plot(tlist, completion_p_displace, color='red', linestyle=symbols[i], linewidth=2, label=r'$\mu$ = '+f'{muList[i]:.2f}')
    plt.plot(tlist, completion_p_squeezing, color='blue', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rList[i]:.2f}')
    
plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title('Passive State', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.show()
    
    
    
    
    
    
'''   
    tlist, Vw_squeezing, completion_squeezing = np.loadtxt(f'./ThermalKinematics/s{r}-cooling.txt', unpack=True, usecols=(0,2,4))
    tlist, Vw_displace, completion_displace = np.loadtxt(f'./ThermalKinematics/d{mu}-heating.txt', unpack=True, usecols=(0,2,4))
    
    

    plt.subplot(121)
    
    plt.plot(tlist, completion_displace, color='red', linestyle=symbols[i], linewidth=2, label=r'$r_d$ = '+f'{muList[i]:.2f}')
    plt.plot(tlist, completion_squeezing, color='blue', linestyle=symbols[i], linewidth=2, label=r'$r_s$ = '+f'{rList[i]:.2f}')
    
    plt.subplot(122)
    
    plt.plot(tlist, completion_p_displace, color='red', linestyle=symbols[i], linewidth=2, label=r'$r_d$ = '+f'{muList[i]:.2f}')
    plt.plot(tlist, completion_p_squeezing, color='blue', linestyle=symbols[i], linewidth=2, label=r'$r_s$ = '+f'{rList[i]:.2f}')


plt.subplot(121)

plt.xscale('log')
plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title('Total State', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.subplot(122)

plt.xscale('log')
plt.xlabel('Time', fontsize=12)
plt.title('Passive State', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()
'''   
    














