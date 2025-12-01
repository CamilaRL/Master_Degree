import numpy as np
import matplotlib.pyplot as plt


def Derivada(xlist, ylist):
    
    ## https://www.youtube.com/watch?v=utRKIlOZbtw
    
    yprime = np.diff(ylist)/np.diff(xlist)
    xprime = []
    
    for i in range(len(yprime)):
        
        xtemp = (xlist[i+1] + xlist[i])/2
        xprime = np.append(xprime, xtemp)
    
    return xprime, yprime


##### MAIN #####

modo = 'Cooling'
dSr_init = 0.1

curvas, cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr_init}/cmod.txt', unpack=True)


cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, len(cmod))))    

dSr = []
dtlist = []
for k in range(len(curvas)):
    
    tlist, Srt = np.loadtxt(f'./FisherInformation_{modo}_{dSr_init}/Sr_curve_{k}.txt', unpack=True)
    
    dt, dSrt = Derivada(tlist, Srt)
    
    for i in range(len(dSrt)):
        
        dSrt[i] = - dSrt[i]
        
    dtlist.append(dt)
    dSr.append(dSrt)
    
    cor = next(colors)
    
    plt.plot(tlist, Srt, color=cor, label=f'|c| = {cmod[k]:.4f}')


plt.legend(loc='center right')
plt.ylabel('Relative Entropy')
plt.xlabel('Time')
plt.title(modo)
plt.show()
    

colors = iter(cmap(np.linspace(0.01, 1, len(cmod))))    

for j in range(len(dtlist)):
    
    cor = next(colors)
    
    plt.plot(dtlist[j], dSr[j], color=cor, label=f'|c| = {cmod[j]:.4f}')
    
plt.ylabel('Entropy Production Rate')
plt.xlabel('Time')
plt.legend()
plt.show()
    