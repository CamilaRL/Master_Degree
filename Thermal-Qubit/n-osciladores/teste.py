import numpy as np
import matplotlib.pyplot as plt
from qutip import *



Tbanho = 10
w = 2
Tqubit = 2
w0 = 2
p = np.exp(w0/(2*Tqubit))/(2*np.cosh(w0/(2*Tqubit)))
gamma = 3

n = 1/(np.exp(w/Tbanho) - 1)

coerencia_max = np.sqrt(max(0.0, p*(1-p)))

tlist = np.arange(0, 10, 0.01)

clist = []
ilist = np.arange(0, 1, 0.1)
for i in ilist:
    for j in ilist:
    
        c = complex(i, j)
        mag = abs(c)
        if mag > coerencia_max:
            
            c = (c/mag) * coerencia_max
        
        clist.append(c)
            
clist = np.sort_complex(clist)


for c in clist:
        
    for t in tlist:
        
        rx = c.real * np.exp(-2*gamma*(n + 0.5)*t)
        ry = -c.imag * np.exp(-2*gamma*(n + 0.5)*t)
        rz = (1/(2*n + 1)) - 2*(p - n/(2*n + 1)) * np.exp(-2*gamma*(2*n + 1)*t)
        
        r_mod = rx**2 + ry**2 + rz**2
        
        coerencia_01 = abs(complex(rx, -ry))
        coerencia_10 = abs(complex(rx, ry))
        
        if coerencia_01 > coerencia_max or coerencia_10 > coerencia_10:
            
            print(f'{c} Não positiva t = {t}')
            print(f'    {coerencia_max} {coerencia_01} {coerencia_10}')