

import numpy as np
from math import *
import random

def trace(M):
	s=0;
	for (i,x) in enumerate(M):
		s=s+M[i,i]
	return s


def form(M):
	M2=matrix(M)
	return lambda A : trace(A*M2.T) 

def formF(n,f):
	M2=matrixF(n,f)
	return lambda A : trace(A*M2.T) 



def Iform(M):
	return lambda A : trace(A*M.T) 


def matrix(a):
	return np.asmatrix(np.array( a));

def matrixF(n,f):
	return matrix( [ [ f(i,j) for j in range(n) ] for i in range(n) ] )



class MatrixLaw:
	def __init__(self,E,M,V,L,size):
		self.size=size
		self.E=E
		self.M=M
		self.M2= matrix( [ [ M[i,j]**2+V[i,j] for j in range(E.shape[1])] for i in range(E.shape[0]) ] )
		self.L=L
		self.Es= [ E**i for i in range(size) ]
		self.Norm= self.L(E*self.Es[size-1])

	def correlation(self,d):
		Er=self.Es[N-d-2]
		Em=self.Es[d]
		res=self.L( self.M* Em *self.M*Er)
		return res/self.Norm
	
	def moy(self,t):
		El=self.Es[t]
		Er=self.Es[self.size-t-1]
		return self.L(El*self.M*Er)/self.L(El*self.E*Er)
	
	def std(self,t):
		El=self.Es[t]
		Er=self.Es[self.size-t-1]
		M2=self.L(El*self.M2*Er)/self.Norm
		M=self.L(El*self.M*Er)/self.Norm
		return sqrt(M2-M*M)


	def correlationT(self,t1,t2):
		if t1==t2:
			return float('nan')
		else:
			tl=min(t1,t2)
			tr=max(t1,t2)
			El=self.Es[tl]
			Em=self.Es[tr-tl-1]
			Er=self.Es[self.size-1-tr]
			norm=self.Norm
			moyL=self.L(El*self.M*Em*self.E*Er)/norm
			moyR=self.L(El*self.E*Em*self.M*Er)/norm
			moy2L=self.L(El*self.M2*Em*self.E*Er)/norm
			moy2R=self.L(El*self.E*Em*self.M2*Er)/norm
			cstd= sqrt( (moy2L-moyL**2)*(moy2R-moyR**2) )
			res=(self.L( El*self.M* Em *self.M*Er))/norm-moyR*moyL

		return res/cstd

	def correlation(self):
		Em=self.Es[0]
		Er=self.Es(self.size-2)
		norm=self.L(Er*self.Es[2])
		R=[]
		for i in range(N):
			R.append(self.L(self.M*Em*self.M*Er)/norm)
			Em=Em*self.E
			Er=self.Es[N-i-2]
		return R


def correlationArr(ML,t):
	return [ML.correlationT(t,i) for i in range(ML.size)]






