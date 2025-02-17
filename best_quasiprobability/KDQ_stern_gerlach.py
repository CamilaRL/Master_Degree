import numpy as np
import matplotlib.pyplot as plt
import cmath

c_real = np.arange(-5, 5, 0.1)
c_img = np.arange(-5, 5, 0.1)

clist = []

for cr in c_real:
	for ci in c_img:
		clist.append(complex(cr, ci))

rho_diag = [0, 1]

def qzx(z, x, rho_diag, rho_c):

	quasi = 0
	
	if z==1 and x ==1:
		quasi = 0.5*(rho_c + rho_diag[0])
		
	if z==-1 and x ==1:
		quasi = 0.5*(np.conj(rho_c) + rho_diag[1])
		
	if z==1 and x ==-1:
		quasi = 0.5*(rho_diag[0] - rho_c)
		
	if z==-1 and x ==-1:
		quasi = 0.5*(rho_diag[1] - np.conj(rho_c))
		
	return quasi.real
	
q11 = []
q_11 = []
q1_1 = []
q_1_1 = []

for rhoc in c_real:

	q11.append(qzx(1, 1, rho_diag, rhoc))
	q_11.append(qzx(-1, 1, rho_diag, rhoc))
	q1_1.append(qzx(1, -1, rho_diag, rhoc))
	q_1_1.append(qzx(-1, -1, rho_diag, rhoc))
	
plt.plot(c_real, q11, label='(1,1)')
plt.plot(c_real, q_11, label='(-1,1)')
plt.plot(c_real, q1_1, label='(1,-1)')
plt.plot(c_real, q_1_1, label='(-1,-1)')
plt.ylabel("Kirkwood-Dirac Quasiprobability")
plt.xlabel("Real Coherent Term")
plt.legend()
plt.show()
