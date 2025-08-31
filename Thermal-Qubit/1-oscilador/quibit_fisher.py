import numpy as np
import matplotlib.pyplot as plt
import qutip as qt


def KDQ_ij(omega, delta, Delta, t):
	
	
	
	qmm_r = (p*(delta**2 + 2*omega**2) - c*delta*omega + delta*(p*delta + c*omega)*np.cos(omega*t))/(2*(Delta**2))
	qmM_r = delta*(p*delta + c*omega)*(np.sin(omega*t/2)*np.sin(omega*t/2))/(Delta**2)
	qMm_r = delta*((1-p)*delta - c*omega)*(1 - np.cos(omega*t))/(2*(Delta**2))
	qMM_r = ((1-p)*(delta**2 + 2*omega**2) + c*delta*omega + delta*((1-p)*delta - c*omega)*np.cos(omega*t))/(2*(Delta**2))
	
	qmm_i = (c*delta*np.sin(omega*t))/(2*Delta)
	qmM_i = -(delta*c*np.cos(omega*t/2)*np.sin(omega*t/2))/Delta
	qMm_i = -(delta*c*np.sin(omega*t))/(2*Delta)
	qMM_i = (c*delta*np.sin(omega*t))/(2*Delta)
	
	qmm = complex(qmm_r, qmm_i)
	qmM = complex(qmM_r, qmM_i)
	qMm = complex(qMm_r, qMm_i)
	qMM = complex(qMM_r, qMM_i)
	
	return qmm, qmM, qMm, qMM
	
	
def KDQ(taulist, omega, delta, Delta):
	
	qmm_list = []
	qmM_list = []
	qMm_list = []
	qMM_list = []
	
	qmm_r_list = []
	qmM_r_list = []
	qMm_r_list = []
	qMM_r_list = []
	
	for tau in taulist:
	
		qmm, qmM, qMm, qMM = KDQ_ij(omega, delta, Delta, tau*np.pi/omega)
		
		qmm_list.append(qmm)
		qmM_list.append(qmM)
		qMm_list.append(qMm)
		qMM_list.append(qMM)
		
		qmm_r_list.append(qmm.real)
		qmM_r_list.append(qmM.real)
		qMm_r_list.append(qMm.real)
		qMM_r_list.append(qMM.real)

	plt.plot(taulist, qmm_r_list, label=r'$q_{--}$')
	plt.plot(taulist, qmM_r_list, label=r'$q_{-+}$')
	plt.plot(taulist, qMm_r_list, label=r'$q_{+-}$')
	plt.plot(taulist, qMM_r_list, label=r'$q_{++}$')

	plt.hlines(y=0, xmin=0, xmax=2, color='black')
	plt.ylabel('Real KDQ')
	plt.xlabel(r'$\omega t / \pi$')
	plt.legend()
	plt.show()
	
	return [qmm_list, qmM_list, qMm_list, qMM_list]


def Extended_KDQ(qKDQ):

	eKDQ_mm = []
	eKDQ_mM = []
	eKDQ_Mm = []
	eKDQ_MM = []
	
	for t in range(len(qKDQ[0])):
	
		sumKDQt = 0
		
		for i in range(4):
			
			sumKDQt = sumKDQt + qKDQ[i][t]
			
		eKDQ_mm.append(qKDQ[0][t]/sumKDQt)
		eKDQ_mM.append(qKDQ[1][t]/sumKDQt)
		eKDQ_Mm.append(qKDQ[2][t]/sumKDQt)
		eKDQ_MM.append(qKDQ[3][t]/sumKDQt)
		
	
	return [eKDQ_mm, eKDQ_mM, eKDQ_Mm, eKDQ_MM]
	
	
def Fisher_Info(eKDQ, taulist, a):
	
	IQt = []
	IQt_r = []
	IQt_i = []
	
	for t in range(len(taulist)):
	
		termo1 = a[0]*a[0]*eKDQ[0][t] + a[0]*a[1]*eKDQ[1][t] + a[1]*a[0]*eKDQ[2][t] + a[1]*a[1]*eKDQ[3][t]
		
		termo2 = (abs( a[0]*(eKDQ[0][t] + eKDQ[1][t]) + a[1]*(eKDQ[2][t] + eKDQ[3][t]) ))**2 
		
		IQt.append(4*(termo1 - termo2))
		IQt_r.append(IQt[t].real)
		IQt_i.append(IQt[t].imag)
		
	
	plt.plot(taulist, IQt_r, label='Real')
	plt.plot(taulist, IQt_i, label='Imaginary')
	plt.ylabel('Quantum Fisher Information')
	plt.xlabel(r'$\omega t / \pi$')
	plt.legend()
	plt.show()

	return IQt


def A_Func(delta, omega):

	A = (delta*qt.sigmaz() + omega*qt.sigmax())/2
	
	evalues = A.eigenenergies()

	return evalues

##### MAIN #####


p = 1/2
c = 1/2
delta = 1
omega = (2**(1/2) - 1)*delta
Delta = (delta**2 + omega**2)**(1/2)

taulist = np.arange(0, 2.1, 0.1)

qKDQ = KDQ(taulist, omega, delta, Delta)

eKDQ = Extended_KDQ(qKDQ)

a = A_Func(delta, omega)

IQt = Fisher_Info(eKDQ, taulist, a)



