
from random import randrange,shuffle,random

def skelSparse(d , p ):
	''' Generate a random digraph with d nodes, p edges and no strongly connected component.  
	The rejection algoritm used slows down considerably if p comes close to the upper limit p <= d*(d-1)/2 '''
	descendants= [frozenset([i]) for i in range(p) ]
	edges = set([]) 
	choice = lambda : randrange (d) 
	while len(edges)< p:
		s,e=  choice(), choice()
		if not( s in desc[e] or (s,e) in edges ) :
			edges.add((s,e))
			descendants[s]=  descendants[s].union(desc[e])
	return edges

def skelFull(d , p ):
	''' Generate a random digraph with d nodes, p edges and no strongly connected component.
	The algorithm used is more memory expansive than the sparse version but its execution time is determinist. '''
	l= range(d)
	full= set(l)
	transitions = [ { j for j in range(i) } for i in range(d)]
	
	def choice(transitions):
		nEdges=sum( ( len(ts) for ts in transitions) ) 
		s= randrange(nEdges) 
		for (i,ts) in enumerate(transitions):
			if len(ts)>s:
				ends= [edge for edge in ts]
				return (i,ends[s])
			else:
				s-=len(ts)

	def chooseEdge():
		chosen= choice(transitions)
		(s,e)= chosen  
		transitions[s].remove(e)
		return chosen
	return [chooseEdge() for i in range (p) ]



def sparseScGraph(d):
	''' Generate a set of edges creating a strongly connected grah 
	    The set is constructed by adding uniformly random edges to the digraph until it is strongly connected  '''
	edges=set() #the set sets of edges
	classesId= [  i for i in range (d) ] # classesId[k] is the class containing the vertex k
	classes= { i:frozenset({})  for i in range (d) } # dictionnary of classes
	descendants= { i:frozenset({}) for i in range(d) } # descendants[k] is the set of all descendants of the class k
	choice= lambda : (randrange(d), randrange(d)) # random selection of a new edge
	while len(classes)>1:
		(s,e)=choice() # selecting new edges
		if not( (s,e) in edges):
			edges.add((s,e))
			sC,eC=classesId[s],classesId[e]
			descendants[sC]=descendants[sC] | descendants[eC] | {e} #updating descendants lists 
			if s in descendants[eC] and not (sC==eC): 
				''' if s was a descendent of the class(e), the classes c(e) and c(s) must fusionned ''' 
				classes[sC] = classes[sC] | classes[eC] | {e,s}
				''' The descendants must be fusionned to '''
				descendants[sC]= descendants[sC] | descendants [eC] | {e,s} 
				''' The we prepare the suppression of class(sd) by updating the index classesId '''
				for n in classes[eC]:
					classesId[n]= classesId[s]					
				classesId[e]=classesId[s]
				''' There should not be any mention of eC left, we can delete the corresponding classes '''
				del classes[eC]
				del descendants[eC]
	''' There is only one strongly connected class at this point '''
	return edges


from sympy.matrices import *
def stoMat(d, dG, L,  gen ):
	''' Generate a random stochastic d×d matrix S with non-zero coefficients along the digraph dG and return c* L 
	The value of the coefficient are chosen accordingly to gen '''
	m=zeros(d,d)
	if (d==1):
		m[0,0]=L
		return m
	s= [ 0 for i in range(d) ]
	for (i,j) in dG:
		m[i,j] =  gen()
		s[i]+= m[i,j]
	for i in range(d):
		for j in range(d):
			m[i,j] = L*m[i,j]/s[i]
	return m 


def randomEdgesT(n,m,p):
	""" Generate p random edges inside a n*m block """
	return [ (randrange(n), randrange(m) ) for i in range(p) ]

def transMat(n, m , edges, gen):
	t=zeros(n,m)
	choice= lambda : (randrange(n), randrange(m))
	t[choice() ] += gen()
	return t	


import sympy as sp
def rationalE():
	""" Generate a random rational for the structure matrix """
	num = 1+ randrange(3)
	#denom = 1+randrange(1)
	return sp.Rational(num)


def rationalQ():
	""" Generate a random rational for the moment matrix"""
	num = randrange(100)
	denom = randrange(100)
	return sp.Rational(num,denom)




def blit(perm, k,l, dims, cdims, source, dest ):
	""" Transfer the matrix 'source' into the ('k','l')-block of the matrix 'dest' applying the permutation 'perm'. 
	'cdims' is the cumulative size of the blocks 
  	"""
	for i in range( dims[k] ):
		for j in range ( dims[l] ):
			dest[ perm[ i + cdims[k]], perm[ j +cdims[l] ] ] = source[i,j]

def sysMats(fG,fT, dims, lambdas, Egen, nQ, Qgen ):
	""" Generate random structure matrix and moment matrix 
		* dims [int] are the size of the blocks
		* lambdas [Real] are the dominant eigenvalue of the block
		* fG  [size -> nEdges]  is the number of irreversible transition between class in function of the number of block
		* nQ is the number of moment matrix
		* Egen, Qgen are the structure and moment matrix coefficient generator
		* fT [size->size->nEdges] is the number of edges in the transition block in function of the size of the block
	"""
	d=sum(dims)
	perm = [ i for i in range(d)]
	shuffle(perm) 
	nb=len(dims)
	g=skelFull( nb, fG (nb) )
	E=zeros((d,d))
	Qs= [ zeros((d,d)) for i in range(nQ) ]
	
	cdims=  [ sum(dims[0:i]) for i in range(d) ]
#	print("dims : {}, cdims : {} ".format(dims,cdims))
	for b in range(nb):
		db=dims[b]
		lb=lambdas[b]
		diagG= sparseScGraph(db)
		block= stoMat(db, diagG, lb, Egen)	
		blit(perm,b,b,dims,cdims,block,E)
		for q in range(nQ):
			qFactor=Qgen()
			qblock=transMat(db,db,diagG, lambda :  qFactor*Egen() ) 
			blit(perm,b,b,dims,cdims,qblock,Qs[q] ) 
	for (k,l) in  g:
		dk,dl=dims[k], dims[l] 
		edges= randomEdgesT(dk, dl, fT(dk,dl) )
		t=transMat(dk,dl,edges, Egen)
#		print("k:{} ,l:{}".format(k,l))
		blit(perm, k,l,dims,cdims,t,E)
#		print(E)
		for q in range(nQ):
			qblock=transMat(dk,dl,edges, Qgen) 
			blit(perm,k,l,dims,cdims,qblock,Qs[q] ) 
	return (E,Qs) 
	 	 

	
def udims(nb, mx ):
	""" Dimension generators | uniform on [1,mx] 
		* nb  number of blocks
		* mx maximal size 
	"""
	if mx==1:
		return [1]*nb
	else:
		return [1+randrange(mx-1) for i in range(nb) ]
	
def dlambdas(nb, p):
	""" eigenvalue generator with a dirac mass on 1"""
	def val():
		u=random()
		if u<p:
			v=sp.Rational(1)
		else:
			v=sp.Rational(1+randrange(9), 10)
		return v
	return [ val() for i in range(nb) ]

	
def system( nb, fG=lambda n : int(0.25*n*n), sBmx=5, pMax=0.5, fT= lambda k,l : int(0.25 *k*l) , nQ=2 ):
	""" Generate structure matrix and moments 
		*nb: number of blocks
		*sBmx: maximal size
		*pMax: weight of the dirac mass in the eigenvalue distribution
		*pT:number of transition	
	"""
	dims=udims(nb, sBmx)
	lambdas=dlambdas(nb,pMax) 
	return sysMats(fG, fT, dims, lambdas, rationalE, nQ, rationalQ )

def system2( dims, lambdas, fG=lambda n : int(0.25*n*n), pMax=0.5, fT= lambda k,l : int(0.25 *k*l) , nQ=2 ):
	nb=len(dims)	
	return sysMats(fG, fT, dims, lambdas, rationalE, nQ, rationalQ )





