import numpy as np
import qutip as qt

def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = qt.bell_state(b) * qt.bell_state(b).dag()
    
    I = qt.Qobj( qt.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
    

def NDQP_ZZ_AA(rho):

	ndqp = []

	zeig = qt.tensor(qt.sigmaz(), qt.sigmaz()).eigenstates()
	aeig = qt.tensor(qt.sigmax(), qt.sigmax()).eigenstates()
	
	projectors_z = []
	projectors_a = []
	
	for v in range(len(zeig[1])):
		
		projectors_z.append(zeig[1][v] * zeig[1][v].dag())
		projectors_a.append(aeig[1][v] * aeig[1][v].dag())
	
	for i in range(len(zeig[0])):
	
		for j in range(len(zeig[0])):
			
			z_perp = qt.qeye(projectors_z[i].dims[0]) - projectors_z[i]
			
			M = projectors_a[j] * projectors_z[i] * qt.Qobj(rho) * z_perp
			
			q = np.trace(M.full())
	
			ndqp.append(q)
			
			print(f'{zeig[0][i]} {aeig[0][j]} {q}')
	
	return ndqp
	

def Write(ff, ndqp_list, p):

	for ndqp in ndqp_list:
		ff.write(f'{p} {ndqp.real} {ndqp.imag}\n')



## MAIN ##


bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1, 0.05)

b = 2

ff = open(f"./results/NDQP/ndqp_zz_xx_{bell_states[b]}.txt", "w")

for p in pList:
	rho = Werner_Density_Matrix(p, bell_states[b])
	
	ndqp = NDQP_ZZ_AA(rho)

	Write(ff, ndqp, p)

ff.close()









