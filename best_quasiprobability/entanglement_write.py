import numpy as np
import qutip as q
	
def Bell_Diagonal_Density_Matrix(e):

    bell_density_matrices = []
    
    for b in ['00', '01', '10', '11']:
        bell_density_matrices.append( q.bell_state(b) * q.bell_state(b).dag() )
        
    rho = q.Qobj(np.zeros((4,4)), dims=[[2,2],[2,2]])
    
    for i in range(4):
        rho = rho + e[i]*bell_density_matrices[i]
        
    return rho
	


def Formation_Entropy(rho):

    C_rho = q.concurrence(rho)

    x = 0.5*(1 + np.sqrt(1 - C_rho**2))
    
    Ef = 0
    if x != 1:
        Ef = -x*np.log2(x) - (1-x)*np.log2(1-x)
    
    return Ef
	
	
def tetrahedron(a1, a2, a3):

	is_bell_state = False
	
	plano1 = a1 - a2 + a3
	plano2 = a1 + a2 - a3
	plano3 = a1 - a2 - a3
	plano4 = a1 + a2 + a3
	
	if plano1 >= -1 and plano2  >= -1 and plano3 <= 1 and plano4 <= 1:
		is_bell_state = True
	
	return is_bell_state
	

# MAIN #

f = open(f'./results/entropy/entanglement_entropy.txt', 'w')


a = np.arange(-1, 1, 0.05)
for a1 in a:
	for a2 in a:
		for a3 in a:

			if tetrahedron(a1, a2, a3):
				
				M = [[1,1,-1,1], [1,-1,1,1], [1,1,1,-1], [1,-1,-1,-1]]
				
				e = np.dot(M, [1, a1, a2, a3])/4
				
				if np.isclose(sum(e), 1):
					
					rho = Bell_Diagonal_Density_Matrix(e)

					S = Formation_Entropy(rho)

					f.write(f'{a1} {a2} {a3} {S}\n')

f.close()

