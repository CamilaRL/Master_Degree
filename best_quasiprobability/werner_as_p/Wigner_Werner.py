import numpy as np
import qutip as q

def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = q.bell_state(b) * q.bell_state(b).dag()
    
    I = q.Qobj( q.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho
    

def MUB():
	
	d = 4
	
	### base MUB retirada do MUB_andreas_2003.pdf
	## each vector is associated with a line in the striation
	M0 = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]], dtype=complex)
	M1 = 0.5*np.array([[1,1,1,1], [1,1,-1,-1], [1,-1,-1,1], [1,-1,1,-1]], dtype=complex)
	M2 = 0.5*np.array([[1,-1,-1.j,-1.j], [1,-1,1.j,1.j], [1,1,1.j,-1.j], [1,1,-1.j,1.j]], dtype=complex)
	M3 = 0.5*np.array([[1,-1.j,-1.j,-1], [1,-1.j,1.j,1], [1,1.j,1.j,-1], [1,1.j,-1.j,1]], dtype=complex)
	M4 = 0.5*np.array([[1,-1.j,-1,-1.j], [1,-1.j,1,1.j], [1,1.j,-1,1.j], [1,1.j,1,-1.j]], dtype=complex)
	
	MUBs = [M0, M1, M2, M3, M4] # each Mi is a different striation
	
	projector = []
	
	for M in MUBs:
	
		Pk = []
		
		for vector in M:
			Pk.append(np.outer(vector, vector))
	
		projector.append(Pk)
	
		
	return MUBs, projector
	
def Line(ak, bk, q, p):

	return ak*q + bk*p	
	

def Translation(c):

	return 


def Lines_per_Striation(k, ak, bk, qlim, plim):

	

	q_list = np.arange(-qlim, qlim + 0.1, 0.1)
	p_list = np.arange(-plim, plim + 0.1, 0.1)
	
	line_c = []
	
	for j in range(len(q_list)):
		c = Line(ak, bk, q_list[j], p_list[j])
		
		if np.iscloe(c, 0):
			line_c.append(c)
			
	
	

def Quantum_Net(Pkj, q, p):

	Q = []
	
	
	
	for k in range(len(ak_list)):
		
		c = Line(a[k], b[k], q, p)
		
	
	

#def probabilities(striations, rho:

	
	

## MAIN ##

ak_list = [0, 1, 1, 1, 2]
bk_list = [1, 0, 1, 2, 1]

#striations, Pkj = MUB()

k = 0
Lines_per_Striation(k, ak_list[k], bk_list[k], 1, 1)



