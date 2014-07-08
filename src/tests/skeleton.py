
import sympy as sym
from sympy.matrices import *
import matplotlib.pyplot as plt

sym.init_printing(use_unicode=True)


d=2
E0 = Matrix( [ [ 1,3] ,  [4,5] ] )
Qs[0] = Matrix( [ [ 8 ,7 ] ,  [9 , 10] ] )
Qs[1] = Matrix( [ [ 1 ,7 ] ,  [5 , 10] ] )

for i in range(d):
	for j in range(d):
		Qs[1][i,j]+= Qs[0][i,j]

import random
import math
def Gs(i,j):
	return lambda : random.gauss(Qs[0][i,j], math.sqrt( Qs[1][i,j]) ) 
import Pymatr.synthesis as Syn


(E0, Qs)  = generation.system(12, 20, sBmx=2)

d=E0.shape[0]
A0 = ones((d,d))

# Jordan-Frobenius decomposition
from Pymatr.connectivity import *

Conn=computeConn(E0)
comps=relabelling(E0)
dims = blockDims(labels)

[E,A,Q1,Q2]=[ relabel(comps, M) for M in [E0,A0,Qs[0],Qs[1] ] ]
Bs = cloneBlocks(dims, E )


# Eigentriples computation
from Pymatr.eigentriples import *
eigens = eigentriples(Bs)
dEigen= dominantEigenvalue(eigens)
dominants =findDominants(eigens)

# Reduced statistics computation
from Pymatr.reduction import *
shadow=shadowTransition(E,dims, dominants)
lshadow=shadowLimit(shadow)
Er= reducedE(dominants,dims,eigens,lshadow)
Ar=reducedA(A,dominants,dims, eigens)
[Q1r,Q2r] = [ reducedStatistic(M,dominants, Bs, dims, eigens) for M in [Q1,Q2] ]


# Limits laws computations
from Pymatr.limits import *
paths=maxPaths(Er)
probP = probPaths(Ar,Er,paths)

# LLN limit distribution
llnLaw=llnLimit(Q1r,paths,probP)


#sigma2s=np.array([0.99,0.01,0])

# CLT computation
xtcl= np.arange(-2,2,0.025)
tclLaw=tclLimit(sigma2s,paths, probP )
#ptcl = tclLaw(xtcl) 


print(" \
matrix E : \n {}, \n \
connection matrix Conn :\n {} \n \
weakly / strongly connected component: \n {} \n \
block dimensions: \n {} \n  \
relabelled matrix E : \n {} \n  \
blocks B(i) : \n {} \n \
list of eigenvalue, left and right eigenvectors : \n {} \n \
dominants blocks : \n {} \n \
shadow transition matrix X : \n {} \n \
limit shadow transition matrix X : \n {} \n \
\n {} \n \
reduced E : \n {} \n \
reduced A : \n {} \n \
reduced Q(1) : \n {} \n \
reduced Q(2) : \n {} \n \
maximal paths : \n {} \n \
path probabilities : \n{}\n \
variances : {} \n \
Law of large number limit distribution: \n {} \n \
  ".format(E0,Conn,labels,dims,E,Bs,eigens,dominants,shadow,lshadow,"="*79,Er,Ar,Q1r,Q2r,paths,probP,sigma2s,llnLaw.latex("s")))

#plt.plot(xtcl,ptcl)
#plt.show()
