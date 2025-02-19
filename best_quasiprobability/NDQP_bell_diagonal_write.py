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

def Bell_Diagonal_Density_Matrix_4(e):
	
	qubit_base = np.array([np.array([1,0]), np.array([0,1])])

	psi1 = (Tensor_Product_1D(qubit_base[0], qubit_base[0]) + Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2)
	psi2 = (Tensor_Product_1D(qubit_base[0], qubit_base[0]) - Tensor_Product_1D(qubit_base[1], qubit_base[1]))/np.sqrt(2)
	psi3 = (Tensor_Product_1D(qubit_base[0], qubit_base[1]) + Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2)
	psi4 = (Tensor_Product_1D(qubit_base[0], qubit_base[1]) - Tensor_Product_1D(qubit_base[1], qubit_base[0]))/np.sqrt(2)
		
	rho = e[0]*np.outer(psi1, psi1) + e[1]*np.outer(psi2, psi2) + e[2]*np.outer(psi3, psi3) + e[3]*np.outer(psi4, psi4)
	
	return rho
	
	
def NDQP_ZZ_AA(rho):

	identity = np.array([[1, 0], [0, 1]])
	pauli_x = np.array([[0, 1],[1, 0]])
	pauli_y = np.array([[0, -1j],[1j, 0]])
	pauli_z = np.array([[1, 0],[0, -1]])
	
	II = Tensor_Product_2D(identity, identity)
	ZZ = Tensor_Product_2D(pauli_z, pauli_z)
	AA = Tensor_Product_2D(pauli_y, pauli_y)
	
	ZZ_perp = II - ZZ
	
	return np.trace(np.dot(AA, np.dot(ZZ, np.dot(rho, ZZ_perp))))

def tetrahedron(a1, a2, a3):

	is_bell_state = False
	
	plano1 = a1 - a2 + a3
	plano2 = a1 + a2 - a3
	plano3 = a1 - a2 - a3
	plano4 = a1 + a2 + a3
	
	if plano1 >= -1 and plano2  >= -1 and plano3 <= 1 and plano4 <= 1:
		is_bell_state = True
	
	return is_bell_state

## MAIN ##

ff = open("./results/NDQP/ndqp_zz_yy_negative.txt", "w")

a = np.arange(-1, 1, 0.05)

for a1 in a:
	for a2 in a:
		for a3 in a:
			
			if tetrahedron(a1, a2, a3):
				
				M = [[1,1,-1,1], [1,-1,1,1], [1,1,1,-1], [1,-1,-1,-1]]
				
				e = np.dot(M, [1, a1, a2, a3])/4
	
				if sum(e) == 1:
				
					rho = Bell_Diagonal_Density_Matrix_4([1, a1, a2, a3])
				
					ndqp_ff = NDQP_ZZ_AA(rho)
				
					if ndqp_ff.real < 0:
						ff.write(f'{a1} {a2} {a3} {ndqp_ff.real} {ndqp_ff.imag}\n')


ff.close()









