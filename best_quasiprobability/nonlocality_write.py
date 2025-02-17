import numpy as np


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
	
	
def Hilbert_Schmidt_Distance(a1, a2, a3):
	# formula retirada de nonlocality-bell-diagonal-zanfardino-2023
	return (1/4)*((a1 - a2)**2 + (a1**2 + a2**2)/2 + (a2 - a2)*(np.sqrt(a1**2 + 4) - np.sqrt(a2**2 + 4)) + 2)
	

# MAIN #

f = open('./results/nonlocality.txt', 'w')


a = np.arange(-1, 1, 0.1)

for a1 in a:
	for a2 in a:
		
			#rho = Bell_Diagonal_Density_Matrix_4([1,a1, a2, a3])
		HS = Hilbert_Schmidt_Distance(a1, a2, 0)
			
		f.write(f'{a1} {a2} {HS}\n')
			
			
f.close()
