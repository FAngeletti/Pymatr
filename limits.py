

import sympy as sym
from sympy.matrices import *


def dPaths(E):
	''' Compute all the directed path appearing in matrix E '''
	d=E.shape[0]
	l=[]
	curr=[]
	def rec(p):
		curr.append(p)
		l.append(curr.copy())
		for i in range(p+1,d):
			if E[p,i]>0:
				rec(i)
		curr.pop()
	for k in range(d): 
		rec(k)
	return l
		
def maxPaths(E):
	''' Find the directed path with maxima lenght '''
	paths=dPaths(E)	
#	print("Path list:{}".format(paths))
	l= max( (len(p) for p in paths) )
	mpaths= [ p for p in paths if (len(p)==l) ]
	return mpaths	

def probPath(A,E,path):
	''' Compute the non-normalized probability of a given path '''
	p=1
	for (s,n) in zip(path, path[1:]) :
		p*= E[s,n]
	p*= A[path[0], path[-1]]
	return p

def probPaths(A,E, paths):
	''' Compute and normalized the probability for all path '''
	ps=[ probPath(A,E, p) for p in paths ]
	tot= sum(ps)
	ps= [ p/tot for p in ps]
	return ps


# LLN limits
from Algo.polytope import VolIntersectionSimplex 
def llnPathLimit(PathMus):
	''' Limit laws for the law of large number on a given path (exact )  '''
	mu0= PathMus[0]
	if len(PathMus)==1:
		return byPieces([mu0,mu0], [0,1,0] )
	normal = Matrix([ pos - mu0 for pos in PathMus[1:] ] )
#	print("Normal : {}, mu0 : {}, PathMus {}".format(normal, mu0, PathMus) )
	offset = - mu0
	return VolIntersectionSimplex(normal,offset) # let's not talk about the magic black box here

from Algo.byPieces import *
def llnLimit(mus,paths, Ppaths):
	''' Limit laws for the law of large number on all paths (exact )  '''
#	print ("Paths: {}".format(paths))
	def musAlongPath(path):
		return [ mus[i] for i in path]
	llaws= [ llnPathLimit(musAlongPath(path)) for path in paths]
	llaw= Ppaths[0] * llaws[0]
	for (p, law ) in zip(Ppaths[1:],llaws[1:]):
		llaw+=p*law
	return llaw





# CLT limits
from Algo.utils import *
def tclPathLimit(sigma2s):
	''' Limit law for the central theorem on a given path associated to variance sigma2s (numerical )'''
	def law(x):
		fun=gaussian(x,sigma2s)
		integral=MarkovianIntegral(simplexAlea(len(sigma2s)), fun) 
		return integral
	return law


def tclLimit(sigma2s,paths, Ppaths):
	''' Limits laws on all paths'''
	def sigma2sAlong(path):
		return [sigma2s[i] for i in path ]
	
	plaws= [ tclPathLimit(sigma2sAlong(path)) for path in paths]
	def law(x):
		return sum( (w*plaw(x) for (plaw, w) in zip(plaws,Ppaths) ) ) 
	return law



