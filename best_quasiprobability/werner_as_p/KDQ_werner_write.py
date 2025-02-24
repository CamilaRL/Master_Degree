import numpy as np
import qutip as q

def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = q.bell_state(b) * q.bell_state(b).dag()
    
    I = q.Qobj( q.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
	

def KDQ_ZZ_AA(rho):

	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	ZZ = np.kron(pauli_z, pauli_z)
	AA = np.kron(pauli_x, pauli_x)
	
	
	return np.trace(np.dot(AA, np.dot(ZZ, rho)))
	
	
	
def KDQ_ZI_AI(rho):

	identity = np.array([[1, 0], [0, 1]])
	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	ZI = np.kron(pauli_z, identity)
	AI = np.kron(pauli_x, identity)
	
	return np.trace(np.dot(AI, np.dot(ZI, rho)))
	
	
	
	
## MAIN ##

bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1, 0.05)

b = 3

ff = open(f"./results/KDQ/kdq_zz_xx_{bell_states[b]}.txt", 'w')
fi = open(f"./results/KDQ/kdq_zi_xi_{bell_states[b]}.txt", 'w')

for p in pList:

	rho = Werner_Density_Matrix(p, bell_states[b])
	
	kdq_aa = KDQ_ZZ_AA(rho)
	
	kdq_ai = KDQ_ZI_AI(rho)

	ff.write(f'{p} {kdq_aa.real} {kdq_aa.imag}\n')
	
	fi.write(f'{p} {kdq_ai.real} {kdq_ai.imag}\n')




ff.close()
fi.close()


