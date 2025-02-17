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

def Bell_Diagonal_Density_Matrix(e1, e2, e3, e4):
			
	qubit_base = np.array([np.array([1,0]), np.array([0,1])])

	psi1 = (Tensor_Product_1D(qubit_base[0], qubit_base[0]) + Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2)
	psi2 = (Tensor_Product_1D(qubit_base[0], qubit_base[0]) - Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2)
	psi3 = (Tensor_Product_1D(qubit_base[0], qubit_base[1]) + Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2)
	psi4 = (Tensor_Product_1D(qubit_base[0], qubit_base[1]) - Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2)
	
	rho = e1*np.outer(psi1, psi1) + e2*np.outer(psi2, psi2) + e3*np.outer(psi3, psi3) + e4*np.outer(psi4, psi4)
	
	return rho


def Bell_Diagonal_Density_Matrix_2(rho0):

	identity = np.array([[1, 0], [0, 1]])
	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	pauli_matrix_list = [identity, pauli_x, pauli_y, pauli_z]
	
	rho = Tensor_Product_2D(identity, identity)
	
	ci = []
	
	for i in range(1, len(pauli_matrix_list), 1):

		sigmai = Tensor_Product_2D(pauli_matrix_list[i], pauli_matrix_list[i])

		ci = np.trace(np.dot(rho0, sigmai))
	
		rho = rho + ci*sigmai
		
	rho = rho/4
	
	return rho
	
	
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


def Bell_Diagonal_Density_Matrix_4(a):

	M = [[1,1,1,-1], [1,1,-1,1], [1,-1,1,1], [1,-1,-1,-1]]

	e = np.dot(M, a)/4
	
	qubit_base = np.array([np.array([1,0]), np.array([0,1])])

	psi1 = (Tensor_Product_1D(qubit_base[0], qubit_base[0]) + Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2)
	psi2 = (Tensor_Product_1D(qubit_base[0], qubit_base[0]) - Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2)
	psi3 = (Tensor_Product_1D(qubit_base[0], qubit_base[1]) + Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2)
	psi4 = (Tensor_Product_1D(qubit_base[0], qubit_base[1]) - Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2)
	
	rho = e[0]*np.outer(psi1, psi1) + e[1]*np.outer(psi2, psi2) + e[2]*np.outer(psi3, psi3) + e[3]*np.outer(psi4, psi4)
	
	return rho
	
	


## MAIN ##

print(Bell_Diagonal_Density_Matrix(0.25,0.25,0.25,0.25))

rho0 = np.ones((4,4))
for i in range(4):
	rho0[i][i] = 0.25

print(Bell_Diagonal_Density_Matrix_2(rho0))

print(Bell_Diagonal_Density_Matrix_3([4,4,4]))
