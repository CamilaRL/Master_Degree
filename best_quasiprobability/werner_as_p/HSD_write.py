import numpy as np


def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = q.bell_state(b) * q.bell_state(b).dag()
    
    I = q.Qobj( q.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
    
	
def Hilbert_Schmidt_Distance(a1, a2, a3):
	# formula retirada de nonlocality-bell-diagonal-zanfardino-2023
	return (1/4)*((a1 - a2)**2 + (a1**2 + a2**2)/2 + (a2 - a2)*(np.sqrt(a1**2 + 4) - np.sqrt(a2**2 + 4)) + 2)


# MAIN #

f = open('./results/HSD/HSD.txt', 'w')

bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1, 0.05)

for p in pList:
	rho = Werner_Density_Matrix(p, bell_states[0])
	
	S = Formation_Entropy(rho)

	output.write(f'{p} {S}\n')


				
f.close()
