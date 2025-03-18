import numpy as np
import qutip as qt


def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = qt.bell_state(b) * qt.bell_state(b).dag()
    
    I = qt.Qobj( qt.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
	

def KDQ_ZZ_AA(rho):

	kdq = []

	zeig = qt.tensor(qt.sigmaz(), qt.sigmaz()).eigenstates()
	aeig = qt.tensor(qt.sigmax(), qt.sigmax()).eigenstates()
	
	projectors_z = []
	projectors_a = []
	
	for v in range(len(zeig[1])):
		
		projectors_z.append(zeig[1][v] * zeig[1][v].dag())
		projectors_a.append(aeig[1][v] * aeig[1][v].dag())
	
	for i in range(len(zeig[0])):
	
		for j in range(len(zeig[0])):
			
			M = projectors_a[j] * projectors_z[i] * qt.Qobj(rho)
			
			q = np.trace(M.full())
	
			kdq.append(q)
			#print(f'{zeig[0][i]} {aeig[0][j]} {q}')
	
	return kdq
	

def KDQ_ZI_AI(rho):

	kdq = []

	zeig = qt.tensor(qt.sigmaz(), qt.qeye(2)).eigenstates()
	aeig = qt.tensor(qt.sigmax(), qt.qeye(2)).eigenstates()
	
	projectors_z = []
	projectors_a = []
	
	for v in range(len(zeig[1])):
		
		projectors_z.append(zeig[1][v] * zeig[1][v].dag())
		projectors_a.append(aeig[1][v] * aeig[1][v].dag())
	
	for i in range(len(zeig[0])):
	
		for j in range(len(zeig[0])):
			
			M = projectors_a[j] * projectors_z[i] * qt.Qobj(rho)
			
			q = np.trace(M.full())
	
			kdq.append(q)
			#print(f'{zeig[0][i]} {aeig[0][j]} {q}')
	
	return kdq

	
def Write_File(ff, kdq, p):

	for q in kdq:
	
		ff.write(f'{p} {q.real} {q.imag}\n')
	
	
## MAIN ##

bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1, 0.01)

b = 3

#ff = open(f"./results/KDQ/kdq_zz_xx_{bell_states[b]}.txt", 'w')
ff = open(f"./results/KDQ/kdq_zi_xi_{bell_states[b]}.txt", 'w')

for p in pList:

	rho = Werner_Density_Matrix(p, bell_states[b])
	
	#kdq_aa = KDQ_ZZ_AA(rho)
	kdq_aa = KDQ_ZI_AI(rho)

	Write_File(ff, kdq_aa, p)
	

ff.close()

