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



### MAIN ###

modoList = ['Cooling', 'Heating']
dSr = 0.1

colors = ['blue', 'red']

for m, modo in enumerate(modoList):
    curvas = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(0), ndmin=1)
    cmod = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/cmod.txt', unpack=True, usecols=(1), ndmin=1)
    
    for i in range(len(curvas)):
        
        tlist, Srt = np.loadtxt(f'./FisherInformation_{modo}_{dSr}/Sr_curve_{i}.txt', unpack=True)
        
        tlist, dSrt = Derivada(tlist, Srt)
        
        plt.plot(tlist, dSrt, color=colors[m], label=f'{modo} - |c| = {cmod[i]:.1f}')


plt.legend(loc='center right')
plt.ylabel('Relative Entropy Derivative')
plt.xlabel('Time')
plt.title('Entropy Production')
plt.show()
    
