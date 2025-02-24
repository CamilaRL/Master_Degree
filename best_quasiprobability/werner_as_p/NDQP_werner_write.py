import numpy as np
import qutip as q

def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = q.bell_state(b) * q.bell_state(b).dag()
    
    I = q.Qobj( q.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
    

def NDQP_ZZ_AA(rho):

	identity = np.array([[1, 0], [0, 1]])
	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	II = np.kron(identity, identity)
	ZZ = np.kron(pauli_z, pauli_z)
	AA = np.kron(pauli_x, pauli_x)
	
	ZZ_perp = II - ZZ
	
	return np.trace(np.dot(AA, np.dot(ZZ, np.dot(rho, ZZ_perp))))

## MAIN ##


bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1, 0.05)

b = 3

ff = open(f"./results/NDQP/ndqp_zz_xx_{bell_states[b]}.txt", "w")

for p in pList:
	rho = Werner_Density_Matrix(p, bell_states[b])
	
	ndqp = NDQP_ZZ_AA(rho)

	ff.write(f'{p} {ndqp.real} {ndqp.imag}\n')

ff.close()









