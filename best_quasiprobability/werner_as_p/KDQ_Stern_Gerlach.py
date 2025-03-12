import numpy as np
import qutip as qt

def Density_Matrix():
    
    rho = [[1, -1j], [2j, 3]]
        
    return rho
	

def KDQ_Z_X(rho):

	f = open('./KDQ_SG.txt', 'w')
	
	z_base = [[1,0], [0,1]]
	x_base = [[(ket0 + ket1)/np.sqrt(2) for ket0, ket1 in zip(z_base[0], z_base[1])], [(ket0 - ket1)/np.sqrt(2) for ket0, ket1 in zip(z_base[0], z_base[1])]]
	
	z_eigvals = qt.sigmaz().eigenenergies()
	x_eigvals = qt.sigmax().eigenenergies()
	
	kdq = []
	
	for i in range(2):
		
		zi = np.outer(z_base[i], z_base[i])
		
		for j in range(2):
		
			xj = np.outer(x_base[j], x_base[j])
			
			q = np.trace(np.dot(xj, np.dot(zi, rho)))
			
			f.write(f'{z_eigvals[i]} {x_eigvals[j]} {q.real} {q.imag}\n')
	
	f.close()
	
	
def KDQ_sigma(rho):

	zeig = qt.sigmaz().eigenstates()
	xeig = qt.sigmax().eigenstates()

	projectors_z = []
	projectors_x = []
	
	for v in range(len(zeig[1])):
		
		projectors_z.append(zeig[1][v] * zeig[1][v].dag())
		projectors_x.append(xeig[1][v] * xeig[1][v].dag())
	
	for i in range(len(zeig[0])):
	
		for j in range(len(zeig[0])):
			
			M = projectors_x[j] * projectors_z[i] * qt.Qobj(rho)
			
			q = np.trace(M.full())
	
			print(f'{zeig[0][i]} {xeig[0][j]} {q}')
	
	
### Main ###

rho = Density_Matrix()

KDQ_sigma(rho)

#KDQ_Z_X(rho)
