import numpy as np
import qutip
import galois as gal


def Map_Continous_To_Discrete(N, plim, qlim):
    
    ## Grid properties
    A = 4*plim*qlim
    side_grid = 2*A/(N*N)
    
    qList = np.arange(-qlim, qlim, 0.1)
    pList = np.arange(-plim, plim, 0.1)
    
    qMapList = []
    pMapList = []
    
    ## map q and p into the grid
    for p in pList:
        for q in qList:
            
            for i in range(N):
    
                p_prime = -plim + i*side_grid
                
                if p_prime <= p and p < (p_prime + side_grid):
                    p_map = i
                
                for j in range(N):
                
                    q_prime = -qlim + j*side_grid
                    
                    if q_prime <= q and q < (q_prime + side_grid):
                        q_map = j
            
            
            if q_map not in qMapList:
                qMapList.append(q_map)
                
            if p_map not in pMapList:
                pMapList.append(p_map)
    
    qMapList = np.sort(qMapList)
    pMapList = np.sort(pMapList)

    return qMapList, pMapList
    

def Lines(N, qlim, plim):
    
    a = [0, 1, 1, 1, 1]
    b = [1, 0, 1, 3, 2]

    qMapList, pMapList = Map_Continous_To_Discrete(N, plim, qlim)
    
    
    ## create lines
    
    striations = []
    
    for k in range(len(a)):
            
        ak = a[k]
        bk = b[k]
        
        c_list = []
        lines_list = []
        print(f'k = {k}')
        for p in pMapList:
            for q in qMapList:
                
                c = ak*q + bk*p
                
                c = c%N
                    
                if c in c_list:
                    lines_list[c_list.index(c)].append([q, p])
                    print(f'{c} {[q,p]}')
                else:
                    c_list.append(c)
                    lines_list.append([[q, p]])
                    print(f'{c} {[q,p]}')
    
        striations.append(lines_list)
    
    return striations


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


def Translation_Coefficients(q_map, p_map):

    ## Associate the grid with the galois field to get the coeficients for translation vector
    GF = gal.GF(2**2)
   
    q = GF(q_map)
    q_coef = q.vector()
    
    p = GF(p_map)
    p_coef = p.vector()
   
    return q_coef, p_coef


def Generalized_Pauli_Matrices(dim):

    phase = np.exp(2j*np.pi/dim)
    
    Z = qutip.Qobj( np.diag([phase**d for d in range(dim)]) )
    
    X = qt.Qobj(np.roll(np.eye(d), -1, axis=1))  # Shift identity matrix
    
    return X, Z
    

def Translation(x, y):

    X, Z = Generalized_Pauli_Matrices(4)
    
    T = (X**(x[0])) * (Z**(y[0]))
    
    for e in range(1, len(x), 1):
    
        T = qutip.tensor(T, (X**x[e])*(Z**y[e]))
        
    
    return T


def Quantum_Net(P, striations):

    Q = []

    for k in range(len(P)):
    
        for j in range(len(P[k])):
            
            xy = striations[k][0]
            
            q_coef, p_coef = Translation_Coefficients(xy[0], xy[1])
            
            T = Translation(q_coef, p_coef)
            
            Q.append( T * P[k][j] * T.dag() )
            
    
    return Q

    
def Werner_Density_Matrix(p, b):
    
    bell_density_matrix = qutip.bell_state(b) * qutip.bell_state(b).dag()
    
    I = qutip.Qobj( qutip.qeye(4)/4 , dims=[[2,2],[2,2]])
    
    rho = p * bell_density_matrix + (1-p) * I
        
    return rho


## MAIN ##

N = 4

plim = 1
qlim = 1


striations = Lines(N, qlim, plim)

#mubs, Pkj = MUB()

#