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

f = open(f'./results/entropy/entropy.txt', 'w')


qubit_base = np.array([np.array([1,0]), np.array([0,1])])

iqb = [ [0,0] , [0,1] , [1,0] , [1,1] ] # Index from Quibit Base

base = []
for i in iqb:
	base.append(np.kron(qubit_base[i[0]], qubit_base[i[1]]))


a = np.arange(-1, 1, 0.1)
for a1 in a:
	for a2 in a:
		for a3 in a:

			if abs(a1) == abs(a2) and abs(a1) == abs(a3):
			
				if tetrahedron(a1, a2, a3):
					
					M = [[1,1,-1,1], [1,-1,1,1], [1,1,1,-1], [1,-1,-1,-1]]
					
					e = np.dot(M, [1, a1, a2, a3])/4
					
					if np.isclose(sum(e), 1):
						
						rho = Bell_Diagonal_Density_Matrix(e)

						S = Formation_Entropy(rho)

						f.write(f'{a1} {a2} {a3} {S}\n')

f.close()

