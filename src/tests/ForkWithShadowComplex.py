


import sympy as sym
from sympy.matrices import *
import matplotlib.pyplot as plt



from Pymatr.generation import *
dims=[2,2,2,2]
nQ=2
half= sym.Rational(1,2)
lambdas=[1,half,1,1]
d=sum(dims)
perm = [ i for i in range(d)]
nb=len(dims)

import random
Egen = lambda : random.randrange(1,5)
Qgen = Egen 


g=[ (0,1), (1,2), (1,3) ] 
E=zeros((d,d))
Qs= [ zeros((d,d)) for i in range(nQ) ]

	
cdims=  [ sum(dims[0:i]) for i in range(d) ]

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
		edges= randomEdgesT(dk, dl, 2 )
		t=transMat(dk,dl,edges, Egen)
#		print("k:{} ,l:{}".format(k,l))
		blit(perm, k,l,dims,cdims,t,E)
#		print(E)
		for q in range(nQ):
			qblock=transMat(dk,dl,edges, Qgen) 
			blit(perm,k,l,dims,cdims,qblock,Qs[q] ) 



print(E)

d=E.shape[0]
A = ones((d,d))

import Pymatr.model as Mod


red= Mod.reduced(A,E, Qs )  

import random
import math
def Gs(i,j):
	return lambda : random.gauss( Qs[0][i,j], 1  )  
import Pymatr.synthesis as Syn

from Pymatr.utils import numerical
L= red.dEigen
nsyn=200
Gen = Syn.MatrixRngOpt(numerical(A),numerical(E/L), Gs, nsyn)

def average():
	s= sum(Gen())
	av= s /nsyn
#	print(" \n sum {},  average: {}\n".format(r, av) )
	return av

lln=red.lln()
import Pymatr.byPieces as Bp
Bp.plot(lln)


import Pymatr.histogram as H
nsample=2000
H.plot(nsample, average)
import matplotlib.pyplot as plt


#plt.xticks( [ Q1r[0].evalf() ], ["µ"] )
plt.show()	



