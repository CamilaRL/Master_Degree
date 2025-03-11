import numpy as np
import qutip as q

def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = q.bell_state(b) * q.bell_state(b).dag()
    
    I = q.Qobj( q.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho


def CHSH(rho):
	
	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	QS = np.kron(pauli_z, (-pauli_z - pauli_x)/np.sqrt(2))
	QT = np.kron(pauli_z, (pauli_z - pauli_x)/np.sqrt(2))
	RS = np.kron(pauli_x, (-pauli_z - pauli_x)/np.sqrt(2))
	RT = np.kron(pauli_x, (pauli_z - pauli_x)/np.sqrt(2))
	
	return np.trace(np.dot(rho, QS)) + np.trace(np.dot(rho, RS)) + np.trace(np.dot(rho, RT)) - np.trace(np.dot(rho, QT))


# MAIN #

bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1.05, 0.05)

b = 0
output = open(f'./results/CHSH/CHSH_{bell_states[b]}.txt', 'w')


for p in pList:
	
	rho = Werner_Density_Matrix(p, bell_states[b])
	
	chsh = CHSH(rho)

	output.write(f'{p} {chsh.real} {chsh.imag}\n')
				
				
output.close()
