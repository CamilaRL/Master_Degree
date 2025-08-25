import numpy as np
import matplotlib.pyplot as plt


def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(-beta*w) - 1)


def Coerencia(c, g, nbar, w, t):
    
    return c*np.exp(4*g*g.conjugate()*(2*nbar + 1)*(np.cos(w*t) - 1)/w**2)
    

### MAIN ###

Tbanho = 2
w = 2
g = complex(0.1, 0.1)

nbar = nbarFunc(Tbanho, w)

curvas = [1 , 1] # aquecer , resfriar
processo = ['Aquecer', 'Resfriar']


for i in range(len(curvas)):

    fig = plt.figure(figsize=(10,5))

    cList = np.loadtxt(f'./FisherInformation_{processo[i]}/c_curva_{curvas[i]}.txt', dtype=complex, unpack=True, ndmin=1)
    time, velocity = np.loadtxt(f'./ThermalKinematics_{processo[i]}/velocity_{curvas[i]}.txt', unpack=True, ndmin=1)
    time, completion = np.loadtxt(f'./ThermalKinematics_{processo[i]}/completion_{curvas[i]}.txt', unpack=True, ndmin=1)

    for c in cList:
        
        coerenciaListReal = []
        coerenciaListImag = []
        
        for t in time:

            coher = Coerencia(c, g, nbar, w, t)

            coerenciaListReal.append(coher.real)
            coerenciaListImag.append(coher.imag)

        plt.subplot(121)
        plt.plot(time, coerenciaListReal, label=f'Coerência Real: {c.real}')
        
        plt.subplot(122)
        plt.plot(time, coerenciaListImag, label=f'Coerência Imaginária: {c.imag}')
        
    
    plt.subplot(121)
    plt.plot(time, velocity, label='Velocidade')
    plt.plot(time, completion, label='Completude')
    plt.xlabel('Time')
    plt.legend()
    
    plt.subplot(122)
    plt.plot(time, velocity, label='Velocidade')
    plt.plot(time, completion, label='Completude')
    plt.xlabel('Time')
    plt.legend()
    
    plt.suptitle(f'{processo[i]} - Curva {curvas[i]}')
    plt.tight_layout()
    plt.show()
    


