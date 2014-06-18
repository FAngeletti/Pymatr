
def blockRange(v):
	''' Start and end indice of each block '''
	w=[ 0 for x in range(len(v)+1)]
	for (i,s) in enumerate(v):
		w[i+1]=w[i]+s
	rangeb=lambda i: (w[i],w[i+1])
	return rangeb

def shadowTransition(E, dims, dominants ):
	''' Shadow transition matrix '''
	S=E.copy()
	rangeb = blockRange(dims)
	for k in dominants:
		(s,e) = rangeb(k)
		for i in range(s,e):
			for j in range(s,e):
				S[i,j]=0
	return S

def shadowLimit(Shadow):
	''' Compute the limit shadow trasition matrix '''
	d=Shadow.shape[0]
	return ( MId(d) - Shadow ).inv()		
		
def subBlocks(M, dims):
	''' Create a function block(k,l) extracting the k-to-l transition block from matrix M according to the block structure defined in dims''' 
	rangeb = blockRange(dims)
	def block(k,l):
		rl,rr = rangeb(k)
		cl,cr = rangeb(l)
		return M[ rl:rr, cl:cr ]
	return block

def reducedE(dominants, dims, eigentriples, lShadow):
	''' Compute the reduced structure matrix '''
	shBlocks=subBlocks(lShadow,dims)
	r=len(dominants)
	Er=zeros(r,r)
	for (i,entry) in zip(range(r),dominants):
		Er[i,i]=1
		for (j,exit) in zip(range(i+1,r),dominants[i+1:]):
			left=eigentriples[entry][2]
			right=eigentriples[exit][1]
			shadow=shBlocks(entry,exit)
			Er[i,j]= left.dot( shadow * right) 
	return Er 


import sympy as sym
from sympy.matrices import *
from Algo.utils import *



def reducedA(A, dominants, dims, eigentriples, lShadow):
	''' Compute the reduced projection matrix '''
	r=len(dominants)
	sht = lShadow.T
	A2 = sht * A * sht
	Ar=zeros(r,r)
	Ablocks=subBlocks(A2, dims)
	for (i,entry) in zip(range(r),dominants):
		for (j,exit) in zip(range(i,r),dominants):
			right=eigentriples[entry][1]
			left=eigentriples[exit][2]
			Ab= Ablocks(entry,exit)
			Ar[i,j]= sum( ewp( Ab, right*left.T) )
	return Ar




def stationaryState( eigentriple ):
	''' Compute the stationary state associated to the eigentriple of a block '''
	left=eigentriple[1]
	right=eigentriple[2]
	norm= left.dot(right) 
	st =  (right* left.T) / norm
	return st

def stationaryTransitions( eigentriple, block ):
	''' stationary transition rate inside a block '''
	L= eigentriple[0]
	left=eigentriple[1]
	right=eigentriple[2]
	norm= left.dot(right)
#	print("block : {} \n, l*r {} \n".format(block,left*right.T ))
	st = ewp( block, (left* right.T)/ (norm * L) )
	return st


def reducedStatistic(P, dominants, blocks, dims, eigentriples):
	''' Compute reduced statistics '''
	r=len(dominants)
	Pr=[0] * r
	subSt= subBlocks(P,dims)
	for (i, pos) in zip(range(r), dominants):
		b=blocks[pos]
		p=stationaryTransitions( eigentriples[pos], b )
		st=subSt(pos,pos) 
#		print( "stationary transition : {} \n statistics: {}\n".format(p, st) )  
		Pr[i]=sum(ewp(p,st))
	return Pr
		
