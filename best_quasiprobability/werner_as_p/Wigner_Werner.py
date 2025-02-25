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
	M0 = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]], dtype=complex)
	M1 = 0.5*np.array([[1,1,1,1], [1,1,-1,-1], [1,-1,-1,1], [1,-1,1,-1]], dtype=complex)
	M2 = 0.5*np.array([[1,-1,-1.j,-1.j], [1,-1,1.j,1.j], [1,1,1.j,-1.j], [1,1,-1.j,1.j]], dtype=complex)
	M3 = 0.5*np.array([[1,-1.j,-1.j,-1], [1,-1.j,1.j,1], [1,1.j,1.j,-1], [1,1.j,-1.j,1]], dtype=complex)
	M4 = 0.5*np.array([[1,-1.j,-1,-1.j], [1,-1.j,1,1.j], [1,1.j,-1,1.j], [1,1.j,1,-1.j]], dtype=complex)
	
	for m in M1:
		for mm in M1:
			print(abs(np.dot(m, mm))**2)
		
	return MUB
	
def striation(MUB):

	#Pk = [ [] for _ in range(N) ]
	
	d = 2
	
	striation = [[]] # list of liens for each striation
	
	for i in range(len(MUB)):
		for s in range(len(MUB)):
			
			kj = np.dot(MUB[i], MUB[s])
			
			if np.isclose(kj**2, 1): ## same striation, same line
			
				for k in range(len(striation)):
					if len(striation[k]) == 0:
						striation[k].append(i)	
			
			if np.isclose(kj**2, 1/d): ## different striation, different line
				
				ki_exist = False
				ks_exist = False
				
				for k in range(len(striation)):
				
					if i in striation[k]:
						ki_exist = True
				
					if s in striation[k]:
						ks_exist = True
				
				if ki_exist == False:
					striation.append([i])
					
				if ks_exist == False:
					striation.append([s])
				
				
			if np.isclose(kj**2, 0): ## same striation, different line
				
				ki_exist = False
				ks_exist = False
				
				for k in range(len(striation)):
				
					if s in striation[k]:
						ks_exist = True
						ks = k
					if i in striation[k]:
						ki_exist = True
						ki = k
				
				if ki_exist and ks_exist: # if i and s exist in same k but are in different k
					
					if ki < ks:
						striation[ki].append(s)
						striation[ks].remove(s)
					elif ks < ki:
						striation[ks].append(i)
						striation[ki].remove(i)
				
				if ki_exist == False:
					striation[ks].append(i)
						
				if ks_exist == False:
					striation[ki].append(s)
			
			print(f'{i} {s} {kj**2:.2f} {striation}')
	
	
	
## MAIN ##

mub = MUB()

for m in mub:
	print(m)
#striation(MUB)
