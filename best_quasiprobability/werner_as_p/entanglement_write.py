import numpy as np
import qutip as q
	
	
def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = q.bell_state(b) * q.bell_state(b).dag()
    
    I = q.Qobj( q.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
	


def Formation_Entropy(rho):

    C_rho = q.concurrence(rho)

    x = 0.5*(1 + np.sqrt(1 - C_rho**2))
    
    Ef = 0
    if x != 1:
        Ef = -x*np.log2(x) - (1-x)*np.log2(1-x)
    
    return Ef
	
	

# MAIN #

bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1, 0.05)

b = 0
output = open(f'./results/entropy/entropy_{bell_states[b]}.txt', 'w')


for p in pList:
	rho = Werner_Density_Matrix(p, bell_states[b])
	
	S = Formation_Entropy(rho)

	output.write(f'{p} {S}\n')



output.close()

