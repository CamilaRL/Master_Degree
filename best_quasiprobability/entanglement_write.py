import numpy as np

def Tensor_Product_1D(vec1, vec2):

	vec = np.zeros(len(vec1)*len(vec2))

	v = 0
	for v1 in range(len(vec1)):
		for v2 in range(len(vec2)):
		
			vec[v] = vec1[v1]*vec2[v2]
		
			v = v + 1
			
	return vec


def Tensor_Product_2D(matrix1, matrix2):

	matrix = np.zeros((len(matrix1)*len(matrix2), len(matrix1)*len(matrix2)), dtype=complex)

	mm = 0
	for i in range(len(matrix1)):
		for j in range(len(matrix1)):
			
			for k in range(len(matrix2)):
				m = 2*i + k
				for l in range(len(matrix2)):
					mm = 2*j + l
					
					matrix[m][mm] = matrix1[i][j]*matrix2[k][l]
				
	return matrix


def Bell_Diagonal_Density_Matrix_3(c):

	identity = np.array([[1, 0], [0, 1]])
	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	pauli_matrix_list = [identity, pauli_x, pauli_y, pauli_z]
	
	rho = Tensor_Product_2D(identity, identity)
	
	for i in range(1, len(pauli_matrix_list), 1):

		sigmai = Tensor_Product_2D(pauli_matrix_list[i], pauli_matrix_list[i])

		rho = rho + c[i-1]*sigmai
		
	rho = rho/4
	
	return rho


def Partial_Trace(rho, base, subsystem_index):
	
	base_A = []
	
	for b in base:
		base_A.append(tuple(b[si] for si in subsystem_index))
		
	base_A  = set(base_A)
	
	base_A_dict = {b:i for i, b in enumerate(base_A)}
	
	rho_A = np.zeros((len(base_A), len(base_A)), dtype=complex)
	
	for i in range(len(base)):
	
		estado_A = tuple(base[i][si] for si in subsystem_index)
		i_estado_A = base_A_dict[estado_A]
		
		for j in range(len(base)):
			
			estado_A_prime = tuple(base[j][sj] for sj in subsystem_index)
			j_estado_A = base_A_dict[estado_A_prime]
			
			system_B_equal = True
			
			for k in range(len(base[i])):
				if k not in subsystem_index:
					if base[i][k] != base[j][k]:
						system_B_equal = False

			if system_B_equal:
				rho_A[i_estado_A][j_estado_A] = rho_A[i_estado_A][j_estado_A] + rho[i][j]
				
	return rho_A


def Entanglement_Entropy(rho_A):

	eigvals = np.linalg.eigvalsh(rho_A)
	eigvals = eigvals[eigvals > 0]
	
	return -np.sum(eigvals * np.log(eigvals))
	
	

# MAIN #

climit = 10
f = open(f'./results/entanglement_entropy{climit}.txt', 'w')

base = []

qubit_base = np.array([np.array([1,0]), np.array([0,1])])

base.append((Tensor_Product_1D(qubit_base[0], qubit_base[0]) + Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2))
base.append((Tensor_Product_1D(qubit_base[0], qubit_base[0]) - Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2))
base.append((Tensor_Product_1D(qubit_base[0], qubit_base[1]) + Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2))
base.append((Tensor_Product_1D(qubit_base[0], qubit_base[1]) - Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2))

c = np.arange(-climit, climit, 0.5)
for c1 in c:
	for c2 in c:
		for c3 in c:

			rho = Bell_Diagonal_Density_Matrix_3([c1,c2,c3])

			rho_A = Partial_Trace(rho, base, [0,1])
			
			S = Entanglement_Entropy(rho_A)

			f.write(f'{c1} {c2} {c3} {S}\n')

f.close()

