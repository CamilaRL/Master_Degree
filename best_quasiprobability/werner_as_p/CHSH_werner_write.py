import numpy as np
import qutip as qt

def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = qt.bell_state(b) * qt.bell_state(b).dag()
    
    I = qt.Qobj( qt.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho


def CHSH(p):
	
	'''a = [1, 0, 0]
	aa = [0, 1, 0]
	b = [1/(np.sqrt(2)), 1/(np.sqrt(2)), 0]
	bb = [1/(np.sqrt(2)), -1/(np.sqrt(2)), 0]
	
	a_sigma = a[0]*qt.sigmax() + a[1]*qt.sigmay() + a[2]*qt.sigmaz()

	b_sigma = b[0]*qt.sigmax() + b[1]*qt.sigmay() + b[2]*qt.sigmaz()
	
	aa_sigma = aa[0]*qt.sigmax() + aa[1]*qt.sigmay() + aa[2]*qt.sigmaz()
	
	bb_sigma = bb[0]*qt.sigmax() + bb[1]*qt.sigmay() + bb[2]*qt.sigmaz()
	
	B = qt.tensor(a_sigma, b_sigma) + qt.tensor(aa_sigma, b_sigma) + qt.tensor(a_sigma, bb_sigma) - qt.tensor(aa_sigma, bb_sigma)

	M_rho = np.trace((qt.Qobj(rho) * B).full())
	
	#M_rho =  (M_rho/2)**2
	
	L = (M_rho - 2)/(2*np.sqrt(2) - 2)'''
	
	L = (np.sqrt((p**2)*3 - p**2) - 1)/(np.sqrt(2) - 1)
	
	L_rho = max([0, L])
	
	return L_rho


# MAIN #

bell_states = ['00', '01', '10', '11']
pList = np.arange(0, 1.01, 0.01)

b = 0
output = open(f'./results/CHSH/CHSH_{bell_states[b]}.txt', 'w')


for p in pList:
	
	rho = Werner_Density_Matrix(p, bell_states[b])
	
	chsh = CHSH(p)

	output.write(f'{p} {chsh.real} {chsh.imag}\n')
				
				
output.close()
