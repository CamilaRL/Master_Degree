import numpy as np
import qutip as q
import galois as gal

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


def Map_Continous_To_Discrete(p, q, N, plim, qlim):
    
    A = 4*plim*qlim
    side_grid = 2*A/(N*N)
    
    for i in range(N):
    
        p_prime = -plim + i*side_grid
        
        if p_prime <= p and p < (p_prime + side_grid):
            p_map = i
        
        for j in range(N):
        
            q_prime = -qlim + j*side_grid
            
            if q_prime <= q and q < (q_prime + side_grid):
                q_map = j
                
    return q_map, p_map


def Line(ak, bk, q, p):

	return ak*q + bk*p	
    

def X(k):
    return k+1
    
def Z(k):
    return (np.exp(2*np.pi*j*k/r)*, dtype=complex)


def Lines_per_Striation(plim, qlim, N):

    GF = gal.GF(2**2)

    striations = []

    ak_list = GF.elements
    bk_list = GF.elements
    
    pList = np.arange(-plim, plim, 0.1)
    qList = np.arange(-qlim, qlim, 0.1)


    qMapList = []
    pMapList = []
    
    ## map q and p into the grid
    for p in pList:
        for q in qList:
            
            q_map, p_map = Map_Continous_To_Discrete(p, q, N, plim, qlim)
            
            if q_map not in qMapList:
                qMapList.append(q_map)
                
            if p_map not in pMapList:
                pMapList.append(p_map)
    
    qMapList = np.sort(qMapList)
    pMapList = np.sort(pMapList)
    

    ## create lines
    for k in range(len(ak_list)):
        
        ## same set [ak, bk] represents one striation k
        ak = ak_list[k]
        bk = bk_list[k]
        
        c_list = []
        lines_list = []
        print(f'k = {k}')
        
        for p in pMapList:
            for q in qMapList:
                
                ## associate a point if the line c
                c = Line(ak, bk, q, p)
                
                c = c%N
                
                ## check if point [q,p] belongs to a existing line c
                if c in c_list:
                    lines_list[c_list.index(c)].append([q, p])
                    print(f'{c} {[q,p]}')
                    
                ## if not, creates another line c
                else:
                    c_list.append(c)
                    lines_list.append([[q, p]])
                    print(f'{c} {[q,p]}')
    
        ## add lines to striation
        striations.append(lines_list)

    return striations


def Quantum_Net(mub, striations):


    for k in range(len(striation)):
    
        dict = mub[k] + striations[k]



## MAIN ##

N = 4

plim = 1
qlim = 1


#mubs, Pkj = MUB()

striations = Lines_per_Striation(plim, qlim, N)




