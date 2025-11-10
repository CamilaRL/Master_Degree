import numpy as np
import matplotlib.pyplot as plt


def nbarFunc(T, w):

	beta = 1/T
	
	return 1/(np.exp(-beta*w) - 1)
	
	
def Coerencia(c, g, nbar, w, t):
    
    return c*np.exp(4*g*g.conjugate()*(2*nbar + 1)*(np.cos(w*t) - 1)/w**2)
    
    
    
### MAIN ###

c = 0.1
Tbanho = 2
w = 2
g = complex(0.1, 0.1)

tlist = np.arange(0, 10, 0.01)

nbar = nbarFunc(Tbanho, w)

coerenciaListReal = []
coerenciaListImag = []

for t in tlist:

    coher = Coerencia(c, g, nbar, w, t)

    coerenciaListReal.append(coher.real)
    coerenciaListImag.append(coher.imag)


fig = plt.figure(figsize=(10,5))

plt.subplot(121)    
plt.scatter(tlist, coerenciaListReal, s=1, label=f'c: {c}')
plt.ylabel('Real Coherence')
plt.xlabel('Time')
plt.legend()

plt.subplot(122)
plt.scatter(tlist, coerenciaListImag, s=1, label=f'c: {c}')
plt.ylabel('Imaginary Coherence')
plt.xlabel('Time')
plt.legend()

plt.tight_layout()
plt.show()


















