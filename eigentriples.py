


import sympy as sym
from sympy.matrices import *
#classical algorithm to compute the eigen pair
def eigenpairN(M_s,eps=1e-12):
	M=sym.N(M_s)
	d=M.shape[0]
	vold=zeros(d,1)
	v=ones(d,1)
	
	while( (v-vold).norm() > eps ) :
	#	print("v:{}, ε:{}".format(v,(v-vold).norm() ) )
		vold=v
		v= M*v
		v= v/ v.norm()
	lambd = v.dot(M*v)/(v.dot (v) )
	return (lambd, v)


def eigenpairS(M):
	''' Compute the rational roots of M and return the greatest one 
	 Assume that the greatest root is rational
	''' 
	vts=M.eigenvects()
	d=M.shape[0]
	basis=0
	l=-1
#	print( "Matrix: {} ".format(M) )
	for (lc, m ,bc) in vts:
		if  sym.im(lc)==0 and ( lc > l ):
#			print(lc)
			basis=bc
			l=lc			
#	print( "eigeinvalue: {} \n eigenvector : {}" .format(l,basis))
	return (l,basis[0])
	


#classical algorithm to compute the eigen pair + left eigeinvector
def eigentriple(M,eps=1e-10):
	''' Compute the dominant eigenvalue and associated right and left eigenvector for matrix with numerical error eps'''
	(l , right) =eigenpairS(M)
	(l, left) = eigenpairS(M.T)
	#print(left)
	n = left.dot(right)
	left=left/n
	return (l,left,right)


def eigentriples(Bs,eps=1e-10):
	''' Compute the dominant eigenvalue and associated right and left eigenvector for each block '''
	return [eigentriple(b,eps) for b in Bs]



def approx(eps):
	''' Admissible numerical error '''
	def sim(x,y):
		return  abs( (x/y) -1 ) < eps 
	return sim

def dominantEigenvalue(eigentriples):
	''' Find the dominant eigenvalue '''
	lambd= max( ( t[0] for t in eigentriples )  )
	return lambd

def findDominants(eigentriples, test=approx(1e-8) ):
	''' Find the indices of the dominant blocks '''
	lambd= dominantEigenvalue(eigentriples)
	dominants=[]
	for (i,t) in enumerate(eigentriples):
		if test(t[0], lambd):
			dominants.append(i)
	return dominants 
